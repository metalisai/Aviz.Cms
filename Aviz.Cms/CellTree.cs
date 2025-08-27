using System.Runtime.InteropServices;
using System.Runtime.CompilerServices;
using MathNet.Numerics.LinearAlgebra;

namespace Aviz.Cms;

using V3 = System.Numerics.Vector3;
using V2 = System.Numerics.Vector2;
using Debug = System.Diagnostics.Debug;

// Cubical marching squares algorithm

//   y
// | 
// | / z
// |/    x
// ------ 

// cube corners and children placement
//  6    7
//   ----   
//  /---/|
// 2|  3|/ 5   (4 not visible)
//  -----      note that the corners are just bit patterns of 0..8
//  0   1         with bit 0 being x, bit 1 y and bit 2 z

// Faces
// 0 - x==0 face
// 1 - x==1 face   face direction = faceId>>1
// 2 - y==0 face      0 - x
// 3 - y==1 face      1 - y
// 4 - z==0 face      2 - z
// 5 - z==1 face

// Face corners and edges
//  2 _____ 3
//   |  2  |       just like cube corners, the corner index
//   |1 0 3|         has bit pattern of coordinates
//  0------- 1      the edges are indexed from bottom clockwise
///   axis1        faces 0,1 - yz     faces 2,3 - zx     faces 4,5 - xy


// Cube edges
//   Cube edges use same indexing order as face corners
//     Thus the 12 edges are grouped into 3 groups of 4. 4 for each direction.
//     (edgeIdx>>2) is edge direction (see direction indexing in Faces)
//   All cube edges always point in positive direction
//     this makes edge subdivision easy and consistent.

// Face children
//
// y
// ---------
// | 2 | 3 |
// |-------|    Same indexing as face corners
// | 0 | 1 |
// --------- x

// Coordinate systems
//   The octree doesn't know anything about global coordinates.
//     Integer coordinates are used to track cell coorinates within octree.
//     These coordinates are indices into hermite data.
//     Edge interections are indexed by dir,dim1,dim2,dim3
//     dim3 always spans the direction of "dir" for sequential memory layout.
//     Cell corners are indexed by x,y,z integers.
//     When new vertices are created, they are created with normalized coordinates within their cell.
//     Global floating point coordinates are resolved after CMS.
//
// The algorithm could've been implemented with less lookup tables.
//    The alternative would be to just duplicate faces and keep track of
//    shared faces. I think the other option would be better now, but initially
//    I started with having edges as reference types as well.

public class CellTree {
    // cosine of angle above which 2D feature is considered flat
    public const float SHARP_FEATURE_ANGLE_THRESHOLD_2D = 0.80f;
    // cosine of angle above which 3D feature is considered flat
    public const float SHARP_FEATURE_ANGLE_THRESHOLD_3D = 0.90f;

    private int size;
    private int rootCell;
    private HermiteData data;

    // arrays of structures of cells, idx is cellId
    private List<Cell> cellPool = new();
    // arrays of structures of cell corner sign bits, idx is cellId
    private List<byte> cellIsInsideBits = new();
    // arrays of structures of cell faces, idx is cellId
    private List<CellFaces> cellFacePool = new();

    // arrays of structures of faces, idx is faceId
    private List<Face> facePool = new();
    // arrays of structures of face edges, idx is faceId
    private List<FaceEdges> faceEdgePool = new();

    // Structure to hold cell context during surface extraction
    struct CellCtx {
        public Cell cell;
        public I3 nodeMin;
        public int nodeSize;
        public HermiteData data;

        public List<Segment> cellSegments = new();
        public HashSet<Vertex> visited = new();
        public HashSet<Vertex> vertices = new();
        public List<int> loop = new();
        public HashSet<Segment> usedSegments = new();

        public CellCtx(HermiteData data) {
            this.data = data;
        }

        public void SetCell(Cell cell, I3 nodeMin, int nodeSize) {
            this.cell = cell;
            this.nodeMin = nodeMin;
            this.nodeSize = nodeSize;
        }
    }

    // Cell normalized segment vertex
    struct LocalVertex {
        public Vertex vertex;
        public V3 cellCoord;
        public V3? normal;
    }

    struct Component {
        public List<LocalVertex> vertices;
        public List<(int i1, int i2, int d)> segments;
        public V3 centralVertex = default;
        public V3 sharpVertex = default;
        public bool isSharp = default;

        public Component() {
            vertices = new();
            segments = new();
        }
    }

    struct SubdivideCtx {
        public HermiteData data;
        public Dictionary<int, (int faceIdx, int cellId, I3 min, int size)> faceSet;
        public List<(int, I3, int)> leafCells;
        public int maxDepth;
    }

    struct ResolvedEdgeVertex {
        public EdgeCoord edge;
        public V3 localPos3;
        public V2? localNormal;
        public V2 localPos2;
    }

    struct MSEdgeVertex {
        public EdgeCoord edge;
        public int edgeIdx; // MS edge index
        public ResolvedEdgeVertex resolved;
    }

    struct MSSegment {
        public MSEdgeVertex v1;
        public MSEdgeVertex v2;
        public V3? sharpFeature; // in cell normalized coordinates
        public V2 sharpFeature2; // on face
    }

    public CellTree(HermiteData data) {
        this.size = data.size;
        this.rootCell = newCell(); 
        InitRoot(rootCell);
        this.data = data;
    }

    public static (V3[] vertices, int[] indices) ExtractSurface(HermiteData data, int initialSubdivisions) {
        var tree = new CellTree(data);
        return tree.ExtractSurface(initialSubdivisions);
    }

    private (V3[] vertices, int[] indices) ExtractSurface(int initialSubdivisions) {
        Subdivide(initialSubdivisions);
        AdaptiveSubdivide(out var faceSet, out var leafCells);

        // segments of all leaf faces
        List<Segment>[] faceSegments = EvaluateFaces(faceSet);
        // storage for components (max 4 per cell)
        Component[] comps = new Component[4];
        for (int i = 0; i < comps.Length; i++) {
            comps[i] = new Component();
        }
        // output vertices and indices
        List<Vertex> cmsVertices = new();
        List<int> indices = new();

        CellCtx cellCtx = new(data);

        foreach (var (cellId, cellMin, cellSize) in leafCells) {
            cellCtx.cellSegments.Clear();
            cellCtx.SetCell(cellPool[cellId], cellMin, cellSize);

            GetCellFaceSegments(cellId, cellCtx.cellSegments, faceSegments);
            var cornerBits = cellIsInsideBits[cellId]; 

            // find and process all components in the cell
            var loops = GetLoops(cellCtx);
            for (int i = 0; i < loops.Length; i++) {
                var loop = loops[i];
                processComponent(loop, in cellCtx, ref comps[i]);
            }

            // check if diagonal case
            //   0 - wont connect components
            //   1 - diagonal case
            //   2 - inverted diagonal case
            int cellCase = 0;
            for (int i = 0; i < 4; i++) {
                byte diag1 = (byte)((0x1 << i) | (0x80 >> i));
                byte diag0 = (byte)~diag1;
                if (cornerBits == diag1) {
                    // case where 2 diagonal corners are inside
                    cellCase = 1;
                } else if (cornerBits == diag0) {
                    // inverted diagonal case
                    cellCase = 2;
                }
            }

            if (cellCase == 0 || !CheckIntersection3D(in comps[0], in comps[1])) {
                // normal triangle fan case
                for (int i = 0; i < loops.Length; i++) {
                    TriangulateFan(in cellCtx, in comps[i], cmsVertices, indices); 
                }
            } else {
                // intersection found
                //   connect components with 'cylinder'
                Debug.Assert(loops.Length == 2);
                TriangulateConnected(in comps[0], in comps[1], cmsVertices, indices);
            }
        }

        // resolve vertex coordinates - compose final mesh data
        var retVertices = new V3[cmsVertices.Count];
        for (int i = 0; i < retVertices.Length; i++) {
            switch (cmsVertices[i]) {
                case EdgeVertex ev:
                    var res = data.GetIntersection(ev.edge);
                    retVertices[i] = res!.Value.p;
                    break;
                case NewVertex nv:
                    retVertices[i] = data.ResolveVertex(nv);
                    break;
            }
        }
        return (retVertices, indices.ToArray());
	}

    private int newCell() {
        int ret = cellPool.Count;
        cellPool.Add(new Cell());
        cellFacePool.Add(new CellFaces());
        cellIsInsideBits.Add(0);
        return ret;
    }

    private int newFace() {
        int ret = facePool.Count;
        facePool.Add(new Face());
        faceEdgePool.Add(new FaceEdges());
        return ret;
    } 

    private void InitRoot(int rootId) {
        ref Cell cell = ref CollectionsMarshal.AsSpan(cellPool)[rootId];
        ref CellFaces cellFaces = ref CollectionsMarshal.AsSpan(cellFacePool)[rootId];
        // initialize edges
        var edges = new EdgeCoord[12];
        for (int i = 0; i < 12; i++) {
            var dir = i >> 2;
            var (es, _) = Lookups.EdgeOffset(i, size);
            edges[i] = new EdgeCoord(es, size, dir);
        }
        // initialize faces
        for (int i = 0; i < 6; i++) {
            var faceId = newFace();
            ref FaceEdges facee = ref CollectionsMarshal.AsSpan(faceEdgePool)[faceId];
            for (int j = 0; j < 4; j++) {
                facee[j] = edges[Lookups.faceEdgeToCubeEdge[i, j]];
            }
            cellFaces[i] = faceId;
        }
    }

    private void AdaptiveSubdivide(out Dictionary<int, (int faceIdx, int cellId, I3 min, int size)> faceSet, out List<(int, I3, int)> leafCells) {
        int maxDepth = (int)Math.Round(Math.Log2(data.size));
        // leaves before subdivision
        List<(int c, I3 m, int s, int d)> leaves = new();
        // final leaf cells
        leafCells = new();
        // leaf faces
        faceSet = new();

        // use struct instead of dragging everything through stack
        SubdivideCtx ctx = new() {
            data = data,
            faceSet = faceSet,
            leafCells = leafCells,
            maxDepth = maxDepth
        };

        void traverse(int cellId, I3 min, int size, int depth) {
            if (!cellPool[cellId].IsLeaf) {
                for (int i = 0; i < 8; i++) {
                    var child = cellPool[cellId].firstChild + i;
                    var co = Lookups.CornerOffsetI(i, size>>1);
                    var childMin = min + co;
                    traverse(child, childMin, size>>1, depth+1);
                }
            } else {
                leaves.Add((cellId, min, size, depth));
            }
        }
        traverse(rootCell, new I3(), size, 0);
        foreach(var (leaf, leafMin, leafSize, depth) in leaves) {
            AdaptiveSubdivide(leaf, leafMin, leafSize, depth, in ctx);
        }
    } 

    private void AdaptiveSubdivide(int cellId, I3 cellMin, int cellSize, int depth, in SubdivideCtx ctx) {
        bool foundAny = false;
        List<V3> allIntr = new ();

        // collect all edge intersections
        for (int i = 0; i < 12; i++) {
            var (es, _) = Lookups.EdgeOffset(i, cellSize);
            var ecoord = new EdgeCoord(cellMin + es, cellSize, i>>2);
            var intrs = data.FindEdgeIntersections(ecoord);
            foreach(var intr in intrs) {
                allIntr.Add(intr.n);
            }
            if (intrs.Length > 1) {
                foundAny = true;
            }
        }

        // get corner sign bits
        byte corners = 0;
        for (int i = 0; i < 8; i++) {
            var offset = Lookups.CornerOffsetI(i, cellSize);
            var cornerOffset = cellMin + offset;
            bool val = data.isInside[cornerOffset.x, cornerOffset.y, cornerOffset.z];
            corners |= (byte)(val ? 1 << i : 0);
        }
        cellIsInsideBits[cellId] = corners;

        // check if needs subdivision
        bool complx = LikelyToContainComplexSurface(allIntr);
        if ((complx && depth < ctx.maxDepth) || (foundAny && depth < ctx.maxDepth)) {
            SubdivideCell(cellId, cellMin, cellSize);
            if (!cellPool[cellId].IsLeaf) {
                for (int i = 0; i < 8; i++) {
                    var child = cellPool[cellId].firstChild + i;
                    var co = Lookups.CornerOffsetI(i, cellSize>>1);
                    var childMin = cellMin + co;
                    AdaptiveSubdivide(child, childMin, cellSize>>1, depth+1, in ctx);
                }
            }
        }
        else
        {
            ctx.leafCells.Add((cellId, cellMin, cellSize));
        }

        // collect leaf faces for later
        //   done here to avoid another traversal
        for (int i = 0; i < 6; i++) {
            var faceId = cellFacePool[cellId][i];
            if (facePool[faceId].IsLeaf) {
                if (!ctx.faceSet.ContainsKey(faceId)) {
                    ctx.faceSet.Add(faceId, (i, cellId, cellMin, cellSize));
                }
            }
        }
    }

    // the heuristic to detect complex surface
    //  checks pairs of normals for large angles
    //  if normals are similar then surface is likely flat
    private static bool LikelyToContainComplexSurface(List<V3> inormals) {
        if (inormals.Count < 2) return false;
        float max = 0.0f;
        for (int i = 0; i < inormals.Count; i++) {
            for (int j = i + 1; j < inormals.Count; j++) {
                var cosAngle = V3.Dot(inormals[i], inormals[j]);
                var dist = 1.0f - cosAngle;
                max = MathF.Max(max, dist);
            }
        }
        // closer to a 0 means more detail
        return max > 0.2f;
    }

    // constant depth subdivision
    //   used to get initial subdivision
    //   otherwise no edges intersect the geometry
    private void Subdivide(int count) {
        void recurse(int cellId, I3 cellMin, int size, int left) {
            if (left == 0) return;
            SubdivideCell(cellId, cellMin, size);
            for (int i = 0; i < 8; i++) {
                var child = cellPool[cellId].firstChild + i;
                var co = Lookups.CornerOffsetI(i, size>>1);
                I3 childMin = cellMin + co;
                recurse(child, childMin, size>>1, left-1);
            }
        }
        recurse(rootCell, new I3(), size, count);
    }

    // run marching squares on all leaf faces
    private List<Segment>[] EvaluateFaces(Dictionary<int, (int faceIdx, int cellId, I3 min, int size)> faceSet) {
        List<(int, (int faceIdx, int cellId, I3 min, int size))> faces = new();
        faces.EnsureCapacity(faceSet.Count);
        foreach (var kvp in faceSet) {
            if (facePool[kvp.Key].IsLeaf) {
                faces.Add((kvp.Key, kvp.Value));
            }
        }

        // Find all unique leaf faces
        // NOTE: leaf nodes in octree are at least 7/8 percent of total nodes, so it's always more efficient to use array instead of hashmap, if ids are dense
		List<Segment>[] faceSegments = new List<Segment>[facePool.Count];

        /*
        // parallel implementation, doesn't seem too efficient
        //   (most time is spent in ExtractSurface anyways)
        System.Collections.Concurrent.ConcurrentBag<List<(int, List<Segment>)>> bag = new();
        var partitioner = System.Collections.Concurrent.Partitioner.Create(0, faces.Count, 1024);
        var options = new System.Threading.Tasks.ParallelOptions();
        options.TaskScheduler = System.Threading.Tasks.TaskScheduler.Default;
        options.MaxDegreeOfParallelism = 6;
        System.Threading.Tasks.Parallel.ForEach(partitioner, options, x => {
            List<(int, List<Segment>)> els = new();
            for (int i = x.Item1; i < x.Item2; i++) {
                List<Segment> segs = new();
                int faceId = faces[i].Item1;
                var (faceIdx, cellId, nodeMin, nodeSize) = faces[i].Item2;
                EvaluateCellFace(cellId, faceEdgePool[faceId], faceIdx, segs, cellIsInsideBits[cellId], nodeMin, nodeSize);
                els.Add((faceId, segs));
            }
            bag.Add(els);
        });
        foreach(var llist in bag) {
            foreach(var (faceId, segs) in llist) {
                faceSegments[faceId] = segs;
            }
        }*/

        // loop through all leaf faces
        foreach (var (faceId, (faceIdx, cellId, nodeMin, nodeSize)) in faces) {
            List<Segment> segs = new();
            EvaluateCellFace(cellId, faceEdgePool[faceId], faceIdx, segs, cellIsInsideBits[cellId], nodeMin, nodeSize);
            faceSegments[faceId] = segs;
        }

		return faceSegments;
    }

    private static void ResolveCellVertex(in CellCtx ctx, ref LocalVertex lv, Vertex cv) {
        if (cv is EdgeVertex mev) {
            var intr = ctx.data.GetIntersectionNormalized(mev.edge, ctx.nodeMin, ctx.nodeSize);
            lv.normal = intr!.Value.n;
            lv.cellCoord = intr.Value.p;
        } else  if(cv is NewVertex mnv) {
            var dif = mnv.cellMin - ctx.nodeMin;
            lv.cellCoord = (mnv.position*mnv.cellSize + dif.ToV3()) / ctx.nodeSize;
        }
    }

    // Triangulate a single component as triangle fan
    private static void TriangulateFan(in CellCtx ctx, in Component comp, List<Vertex> outV, List<int> outI) {
        int start = outV.Count;
        var vspan = CollectionsMarshal.AsSpan(comp.vertices);
        foreach(ref LocalVertex v in vspan) {
            outV.Add(v.vertex);
        }
        int centerIdx = outV.Count;
        var apex = comp.isSharp ? comp.sharpVertex : comp.centralVertex;
        outV.Add(new NewVertex(apex, ctx.nodeMin, ctx.nodeSize));

        foreach(var (s, e, _) in comp.segments) {
            outI.Add(start + s);
            outI.Add(start + e);
            outI.Add(centerIdx);
        }
    }

    // check if two components intersect in 3D
    //   used to resolve internal ambiguity
    private static bool CheckIntersection3D(in Component comp1, in Component comp2) {
        var apex1 = comp1.sharpVertex;
        var apex2 = comp2.sharpVertex;
        // test all triangle pairs
        for (int i = 0; i < comp1.segments.Count; i++) {
            var (s1, e1, _) = comp1.segments[i];
            var v11 = comp1.vertices[s1].cellCoord;
            var v12 = comp1.vertices[e1].cellCoord;
            for (int j = i + 1; j < comp2.segments.Count; j++) {
                var (s2, e2, _) = comp2.segments[j];
                var v21 = comp2.vertices[s2].cellCoord;
                var v22 = comp2.vertices[e2].cellCoord;
                bool intersect = Util.TriTriIntersect(apex1, v11, v12, apex2, v21, v22);
                if (intersect) { 
                    return true;
                }
            }
        }
        return false;
    }

    // NOTE: not sure if this is the best cost function
    //   cost to minimize triangle size
    //   the cost function only improves triangle quality
    private static float TriangleCost(V3 v1, V3 v2, V3 v3) {
        var e1 = v2 - v1;
        var e2 = v3 - v2;
        var e3 = v1 - v3;
        return V3.Dot(e1, e1) + V3.Dot(e2, e2) + V3.Dot(e3, e3);
    } 

    // Triangulate two connected components, the alternate triangulation
    //   for diagonal ambiguity. Connecting the diagonal loops with 'cylinder'
    // Just dynamic programming didn't seem to be enough, so I had to improvise
    // Connecting line segments is divided into 3 runs
    //   only segments from a pair of faces are connected each run
    //   plain DP requires locally optimal decisions, which produced self-intersections
    //   Maybe I'm missing some magic cost function?
    // This '3-run' constraint assumes that the components 'wrap' diagonal corners
    //   which guarantees that the components dont share any face
    //   and contain segments from exactly 3 faces each
    // Only 2 offsets need to be considered
    //   because connecting segments that belong to
    //   opposite faces is obviously wrong (twists).
    // Might also perform better, instead of n^3 its O(a^2+b^2+c^2) where n=a+b+c
    // Downside - I don't know if this guarantees no self-intersections
    //   (but DP didn't seem to either)

    private static void TriangulateConnected(in Component comp1, in Component comp2, List<Vertex> outV, List<int> outI) {
        // take all vertices except the central one
        var cv1 = comp1.vertices;
        var cv2 = comp2.vertices;
        var verts1 = comp1.segments.Select(x => cv1[x.i1].cellCoord).ToArray();
        var verts2R = comp2.segments.Select(x => cv2[x.i2]).Reverse();
        var verts2 = verts2R.Select(x => x.cellCoord).ToArray();
        var planes1 = comp1.segments.Select(x => x.d).ToArray();
        var planes2 = comp2.segments.Select(x => x.d).Reverse().ToArray();

        int[] runs1 = Util.RunLengths(planes1).ToArray();
        int[] runs2 = Util.RunLengths(planes2).ToArray();

        // expecting that segments will be on exactly 3 planes each
        //   otherwise the algorithm needs to be changed
        Debug.Assert(runs1.Length == runs2.Length && runs1.Length == 3);
        int[,] bestOptions = new int[verts1.Length+1, verts2.Length+1];
        List<bool> bestPath = new();
        float bestCost = float.MaxValue;
        int bestOffset = 0;
        int l2Ofst = 0;
        // try all offsets
        for(int cr = 0; cr < 3; cr++) {
            // can't connect to opposite plane without twist
            if (planes1[0] == planes2[l2Ofst%planes2.Length]) {
                l2Ofst += runs2[cr];
                continue;
            }

            // NOTE: intentionally reused between 3 runs
            //   each run uses subset of the state space
            float?[,] memo = new float?[verts1.Length + 1, verts2.Length + 1];

            // DP to optimize triangle quality
            float next(int li, int lj, int l1Len, int l2Len) {
                if (li >= l1Len && lj >= l2Len) {
                    return 0.0f;
                }

                if (memo[li, lj] is float val) {
                    return val;
                }

                float costA = float.MaxValue;
                float costB = float.MaxValue;
                if (li < l1Len) {
                    var v1 = verts1[li%verts1.Length];
                    var v2 = verts1[(li+1)%verts1.Length];
                    var v3 = verts2[(l2Ofst+lj)%verts2.Length];
                    costA = TriangleCost(v1, v2, v3) + next(li+1, lj, l1Len, l2Len);
                }
                if (lj < l2Len) {
                    var v1 = verts1[li%verts1.Length];
                    var v2 = verts2[(l2Ofst+lj)%verts2.Length];
                    var v3 = verts2[(l2Ofst+lj+1)%verts2.Length];
                    costB = TriangleCost(v1, v2, v3) + next(li, lj+1, l1Len, l2Len);
                }
                var ret = float.Min(costA, costB);
                bestOptions[li, lj] = ret == costA ? 0 : 1;
                memo[li, lj] = ret;
                return ret;
            }

            // divide pairs of faces into 3 runs
            //   pair - for example xy plane from loop1 an yz plane from loop2
            //      then repeat for yz and zx, then zx and xy
            //      however the components are wound
            //   not sure if this prevents self-intersections in all cases
            //   but seems better than any cost function I could come up with
            float cost = 0.0f;
            int gi = 0, gj = 0;
            for (int r = 0; r < 3; r++) {
                memo = new float?[verts1.Length + 1, verts2.Length + 1];
                cost += next(gi, gj, gi+runs1[r], gj+runs2[r]);
                gi += runs1[r];
                gj += runs2[r];
            }

            // store decision sequence if best so far
            if (cost < bestCost) {
                bestPath.Clear();
                int mi = 0, mj = 0;
                while (mi < verts1.Length || mj < verts2.Length) {
                    var opt = bestOptions[mi, mj];
                    bestPath.Add(opt == 0);
                    if (opt == 0) {
                        mi++;
                    } else {
                        mj++;
                    }
                }
                bestCost = cost;
                bestOffset = l2Ofst;
            }
            l2Ofst += runs2[cr];
        }

        // add vertices to output
        int comp1Base = outV.Count;
        foreach(var v in comp1.vertices) {
            outV.Add(v.vertex);
        }
        int comp2Base = outV.Count;
        foreach(var v in verts2R) {
            outV.Add(v.vertex);
        }

        // make sure no edges are missed
        Debug.Assert(bestPath.Count(x => x) == verts1.Length);
        Debug.Assert(bestPath.Count(x => !x) == verts2.Length);

        int i = 0, j = 0;
        // reconstruct best known triangulation
        foreach(var option in bestPath) {
            if(option) {
                outI.Add((comp1Base + i)%verts1.Length);
                outI.Add(comp1Base + ((i+1)%verts1.Length));
                outI.Add(comp2Base + (bestOffset + j)%verts2.Length);
                i += 1;
            } else {
                outI.Add((comp1Base + i)%verts1.Length);
                outI.Add(comp2Base + (bestOffset + j + 1)%verts2.Length);
                outI.Add(comp2Base + (bestOffset + j)%verts2.Length);
                j += 1;
            }
        }
    }

    // resolve component vertices and calculate central and sharp vertex
    private static void processComponent(int[] loop, in CellCtx ctx, ref Component outComp)
    {
        outComp.vertices.Clear();
        outComp.segments.Clear();
        V3 massPoint = V3.Zero;
        int vIdx = 0;
        Dictionary<Vertex, int> vertexMap = new();
        float mass = 0.0f;
        foreach (var idx in loop) {
            var seg = ctx.cellSegments[idx];
            int v1Idx;
            if (!vertexMap.TryGetValue(ctx.cellSegments[idx].v1, out v1Idx)) {
                LocalVertex mv1 = default;
                mv1.vertex = ctx.cellSegments[idx].v1;
                ResolveCellVertex(ctx, ref mv1, seg.v1);
                outComp.vertices.Add(mv1);
                massPoint += mv1.cellCoord;
                v1Idx = vIdx++;
                vertexMap[ctx.cellSegments[idx].v1] = v1Idx;
                mass += 1.0f;
            }
            int v2Idx;
            if (!vertexMap.TryGetValue(ctx.cellSegments[idx].v2, out v2Idx)) {
                LocalVertex mv2 = default;
                mv2.vertex = ctx.cellSegments[idx].v2;
                ResolveCellVertex(ctx, ref mv2, seg.v2);
                outComp.vertices.Add(mv2);
                massPoint += mv2.cellCoord;
                v2Idx = vIdx++;
                vertexMap[ctx.cellSegments[idx].v2] = v2Idx;
                mass += 1.0f;
            }
            outComp.segments.Add((v1Idx, v2Idx, seg.nDir));
        }
        massPoint *= 1.0f / mass;
        outComp.centralVertex = massPoint;

        var vcount = outComp.vertices.Count;
        float minCos = 1.0f;
        for (int i = 0; i < vcount; i++) {
            for (int j = i + 1; j < vcount; j++) {
                var n1 = outComp.vertices[i].normal;
                var n2 = outComp.vertices[j].normal;
                if (n1 == null || n2 == null) continue;
                var cos = V3.Dot(n1.Value, n2.Value);
                if (cos < minCos) {
                    minCos = cos;
                }
            }
        }

        if (minCos < SHARP_FEATURE_ANGLE_THRESHOLD_3D) {
            var withNormalCount = outComp.vertices.Count(v => v.normal != null);
            var matrix = Matrix<float>.Build.Dense(withNormalCount, 3);
            var vector = Vector<float>.Build.Dense(withNormalCount, 1.0f);
            int vi = 0;
            for (int i = 0; i < vcount; i++) {
                if (outComp.vertices[i].normal is V3 n) {
                    matrix.SetRow(vi, new float[] { n.X, n.Y, n.Z });
                    vector[vi] = V3.Dot(n, outComp.vertices[i].cellCoord - massPoint);
                    vi++;
                } else {
                    continue;
                }
            }
            var svd = matrix.Svd(true);
            var mw = svd.W;
            var ms = svd.S;
            // regularization tuning parameter
            float lambda = 1e-2f;
            for (int i = 0; i < ms.Count; i++)
            {
                // ridge regression style regularization
                //   helps counter numberical instability
                ms[i] = ms[i] / (ms[i]*ms[i] + lambda);
                // truncation style regularization
                //   seemed less stable than ridge regression
                //ms[i] = MathF.Abs(ms[i]) > 1e-3 ? 1/ms[i] : 0.0f;
            }
            mw.SetDiagonal(ms);
            var pseudoInverse = (svd.U * mw * svd.VT).Transpose();
            var result = pseudoInverse * vector;
            var p = new V3(result[0], result[1], result[2]) + massPoint;
            // NOTE: not sure how this is best handled
            //   this discards sharp features that go far beyond cell
            //     still use sharp features for ambiguity resolution
            //   a better solution would find the nullspace and cell intersection
            //   then find closest point on that to the mass point
            if (p.X < -0.5f || p.X > 1.5f ||
                p.Y < -0.5f || p.Y > 1.5f ||
                p.Z < -0.5f || p.Z > 1.5f) {
                outComp.isSharp = false;
            } else {
                outComp.isSharp = true;
            }
            outComp.sharpVertex = p;
        } else {
            outComp.isSharp = false;
            outComp.sharpVertex = massPoint;
        }
    } 

    // subdivide a cell into 8 children
    private void SubdivideCell(int cellId, I3 min, int size) {
        Debug.Assert(size > 1 && (size&1) == 0);
        // subdivide faces first
        for (int i = 0; i < 6; i++) {
            var faceId = cellFacePool[cellId][i];
            if (facePool[faceId].IsLeaf) {
                SubdivideFace(faceId, min, size, i);
            }
        }

        // create new children
        int firstChild = newCell();
        // NOTE: we always assume all 8 children are allocated continuously
        for (int i = 0; i < 7; i++) {
            int res = newCell();
            Debug.Assert(res == firstChild + i + 1);
        }
        // NOTE: do not move this above newCell() calls!
        ref Cell cell = ref CollectionsMarshal.AsSpan(cellPool)[cellId];
        cell.firstChild = firstChild;

        // create new edges in the middle of the cell and assign them to children
        // assign outer edges created by face subdivision to children
        ConnectChildFaces(cellId, min, size);
    }

    private void ConnectChildFaces(int cellId, I3 min, int size) {
        var cellSpan = CollectionsMarshal.AsSpan(cellPool);
        var cellFaceSpan = CollectionsMarshal.AsSpan(cellFacePool);
        ref Cell cell = ref cellSpan[cellId];
        ref CellFaces cellFaces = ref cellFaceSpan[cellId];
        Debug.Assert(size > 1 && (size&1) == 0);
        int hsize = size >> 1;
        // outer faces 
        for (int i = 0; i < 24; i++) {
            int child = cell.firstChild + Lookups.childOuterFaceMap[i, 0];
            int fidx = Lookups.childOuterFaceMap[i, 1];
            int face = i >> 2;
            cellFaceSpan[child][face] = facePool[cellFaces[face]].firstChild + fidx;
        }

        // create inner faces (without edges assigned)
        var newFaces = new int[12];
        for (int i = 0; i < 12; i++) {
            newFaces[i] = newFace();
        }
        var faceEdgeSpan = CollectionsMarshal.AsSpan(faceEdgePool);

        for (int i = 0; i < 24; i++) {
            var child = Lookups.childInnerFaceMap[i, 0];
            var childFace = Lookups.childInnerFaceMap[i, 1];
            var newFace = i>>1;
            cellFaceSpan[cell.firstChild + child][childFace] = newFaces[newFace];
        }

        // create new edges at the centers of 4 cells at each of the 6 faces
        var cellCenter = min + new I3(hsize, hsize, hsize);
        EdgeCoord[] newEdges = [
            // x
            new EdgeCoord(min + new I3(0, hsize, hsize), hsize, 0),
            new EdgeCoord(cellCenter, hsize, 0),
            // y
            new EdgeCoord(min + new I3(hsize, 0, hsize), hsize, 1),
            new EdgeCoord(cellCenter, hsize, 1),
            // z
            new EdgeCoord(min + new I3(hsize, hsize, 0), hsize, 2),
            new EdgeCoord(cellCenter, hsize, 2),
        ];

        //                            _|/_
        // assign new edges to faces   /|
        for (int i = 0; i < 24; i++) {
            var newFace = i >> 1;
            var edge = Lookups.faceInnerEdgeMap[i, 0];
            var newEdge = Lookups.faceInnerEdgeMap[i, 1];
            faceEdgeSpan[newFaces[newFace]][edge] = newEdges[newEdge];
        }

        EdgeCoord[,] newInternalEdges = new EdgeCoord[6, 4];
        var neSpan = MemoryMarshal.CreateSpan<EdgeCoord>(ref Unsafe.As<byte, EdgeCoord>( ref MemoryMarshal.GetArrayDataReference(newInternalEdges)), newInternalEdges.Length);
        for (int i = 0; i < 6; i++) {
            FaceInternalEdges(i, size, min, neSpan.Slice(i*4));
        }

        for (int i = 0; i < 24; i++) {
            var face = i >> 1;
            var edge = Lookups.faceOuterEdgeMap[i, 0];
            var iface = Lookups.faceOuterEdgeMap[i, 1];
            var iedge = Lookups.faceOuterEdgeMap[i, 2];
            faceEdgeSpan[newFaces[face]][edge] = newInternalEdges[iface, iedge];
        }
    } 

    // check 2D segment intersection
    //   used to resolve face ambiguity
    private static bool CheckIntersection(MSSegment[] segments) {
        Debug.Assert(segments.Length == 2);
        (V2 s, V2 e)[] checkSegs = new (V2, V2)[4];
        int segCount = 0;
        // collect all segments
        foreach(var seg in segments) {
            if (seg.sharpFeature.HasValue) {
                checkSegs[segCount++] = (
                    seg.v1.resolved.localPos2,
                    seg.sharpFeature2
                );
                checkSegs[segCount++] = (
                    seg.sharpFeature2,
                    seg.v2.resolved.localPos2
                );
            } else {
                checkSegs[segCount++] = (
                    seg.v1.resolved.localPos2,
                    seg.v2.resolved.localPos2
                );
            }
        }
        // check all pairs of segments
        for (int i = 0; i < segCount; i++) {
            for (int j = i + 1; j < segCount; j++) {
                if (Util.SegmentsIntersect(checkSegs[i].Item1, checkSegs[i].Item2, checkSegs[j].Item1, checkSegs[j].Item2)) {
                    return true;
                }
            }
        }
        return false;
    }

    // resolve vertex 2D
    private static void ResolveEdgeVertex(HermiteData data, ref MSEdgeVertex ev, int faceId) {
        var intr = data.GetIntersectionWithLeaf(ev.edge);
        ref ResolvedEdgeVertex res = ref ev.resolved;
        res.edge = intr!.Value.leaf;
        res.localNormal = Lookups.CubeDirToFaceDir(faceId, intr.Value.n);
        int cubeEdge = Lookups.faceEdgeToCubeEdge[faceId, ev.edgeIdx];
        var (e1s, e1e) = Lookups.EdgeOffset(cubeEdge, 1);
        res.localPos3 = V3.Lerp(e1s.ToV3(), e1e.ToV3(), intr.Value.t);
        res.localPos2 = Lookups.CubePosToFacePos(faceId, res.localPos3);
    }

    // run marching squares on a leaf face
    private void EvaluateCellFace(int cellId, FaceEdges edges, int faceId, List<Segment> segments, byte cornerBits, I3 min, int size) {	

        int caseId = 0;
		for (int i = 0; i < 4; i++) {
			int corner = Lookups.faceCubeCorners[faceId, i];
            caseId |= (((cornerBits >> corner) & 1)) << i;
		}

		if (caseId == 0 || caseId == 15) {
			return; // no segments to draw
		}

        static void tryFindSharpFeature(ref MSSegment seg, int faceId) {
            ref ResolvedEdgeVertex ev1 = ref seg.v1.resolved;
            ref ResolvedEdgeVertex ev2 = ref seg.v2.resolved;
            if (ev1.localNormal is V2 n1 && ev2.localNormal is V2 n2) {
                // first line equation (perpendicular to normal)
                var t1 = new V2(n1.Y, -n1.X);
                // second line equation (perpendicular to normal)
                var t2 = new V2(-n2.Y, n2.X);
                // determinant for 2D matrix inversion
                var det = t1.X * t2.Y - t1.Y * t2.X;
                if (float.Abs(det) < 1e-4
                    || V2.Dot(n1, n2) >= SHARP_FEATURE_ANGLE_THRESHOLD_2D) {
                    return;
                }
                var p1 = ev1.localPos2;
                var p2 = ev2.localPos2;
				var t = (t2.Y*(p2.X-p1.X) - t2.X*(p2.Y-p1.Y)) / det;
				var p3 = p1 + t * t1;
				//p3.X = float.Clamp(p3.X, 0.0f, 1.0f);
				//p3.Y = float.Clamp(p3.Y, 0.0f, 1.0f);
                // NOTE: not sure if this is best handled
                //   clamp or just ignore sharp feature?
                //   this discards sharp features that go far beyond cell
                if (p3.X > -0.5f && p3.X < 1.5f
                    && p3.Y > -0.5f && p3.Y < 1.5f) {
                    seg.sharpFeature2 = p3;
                    seg.sharpFeature = Lookups.FacePosToCubePos(faceId, p3);
                }
            }
        }

        var segs = new MSSegment[2];
        segs[0] = new MSSegment();
        segs[1] = new MSSegment();

        // detect sharp feature and solve ambiguous cases
        int segCount = 0;
        for (int alt = 0; alt < 2; alt++) {
            segCount = 0;
            // use alt case if ambiguous and first fails
            int curCase = alt == 0 ? caseId : caseId switch {
                6 => 16, // 6 -> 6b
                9 => 17, // 9 -> 9b
                _ => 0
            };
            for (int i = 0; i < 4 && Lookups.quadSegments[curCase, i] >= 0; i+=2) {
                ref MSSegment seg = ref segs[i>>1];

                seg.v1.edgeIdx = Lookups.quadSegments[curCase, i];
                seg.v1.edge = edges[seg.v1.edgeIdx];
                ResolveEdgeVertex(data, ref seg.v1, faceId);

                seg.v2.edgeIdx = Lookups.quadSegments[curCase, i + 1];
                seg.v2.edge = edges[seg.v2.edgeIdx];
                ResolveEdgeVertex(data, ref seg.v2, faceId);
                tryFindSharpFeature(ref seg, faceId);
 
                segCount++;
            }
            // should never try alt case with only 1 segment
            Debug.Assert(alt == 0 || segCount == 2);

            // only 1 segment - not ambiguous
            if (segCount == 1) {
                break;
            } else {
                Debug.Assert(segCount == 2, "Should be 1 or 2 segments");
                // check if segments intersect
                if (CheckIntersection(segs)) {
                    if (alt == 0) {
                        // segments intersect, so we need to use alternative case
                        continue;
                    } else {
                        // somehow segments still intersect
                        //   just discard sharp features
                        segs[0].sharpFeature = null;
                        segs[1].sharpFeature = null;
                    }
                }
            }
        }

        // add segments
        for (int i = 0; i < segCount; i++) {
            ref MSSegment seg = ref segs[i];
            if (!seg.sharpFeature.HasValue) {
                segments.Add(new Segment(
                    new EdgeVertex(seg.v1.resolved.edge),
                    new EdgeVertex(seg.v2.resolved.edge),
                    faceId>>1));
            } else {
                var sfv = seg.sharpFeature.Value;
                var newVert = new NewVertex(sfv, min, size);
                var toNewVert1 = seg.sharpFeature2 - seg.v1.resolved.localPos2;
                var toNewVert2 = seg.sharpFeature2 - seg.v2.resolved.localPos2;
                float minDist = float.Max(V2.Dot(toNewVert1, toNewVert1), V2.Dot(toNewVert2, toNewVert2));
                // if new vertex is very close to one of the vertices
                //   then just discard it
                //   don't want to create very small triangles
                if (minDist < 1e-3) {
                    segments.Add(new Segment(
                        new EdgeVertex(seg.v1.resolved.edge),
                        new EdgeVertex(seg.v2.resolved.edge),
                        faceId>>1));
                } else {
                    segments.Add(new Segment(
                        new EdgeVertex(seg.v1.resolved.edge),
                        newVert,
                        faceId>>1));
                    segments.Add(new Segment(
                        newVert,
                        new EdgeVertex(seg.v2.resolved.edge),
                        faceId>>1));
                }
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void FaceInternalEdges(int faceIndex, int cellSize, I3 cellMin, Span<EdgeCoord> outArray) {
        var hsize = cellSize >> 1;
        var center = faceIndex switch {
            0 => new I3(0, hsize, hsize),
            1 => new I3(cellSize, hsize, hsize),
            2 => new I3(hsize, 0, hsize),
            3 => new I3(hsize, cellSize, hsize),
            4 => new I3(hsize, hsize, 0),
            5 or _ => new I3(hsize, hsize, cellSize),
        };
        center = cellMin + center;
        var dir = faceIndex>>1;
        var axis1 = dir switch {
            0 => new I3(0, hsize, 0),
            1 => new I3(0, 0, hsize),
            2 or _ => new I3(hsize, 0, 0),
        };
        var axis2 = dir switch {
            0 => new I3(0, 0, hsize),
            1 => new I3(hsize, 0, 0),
            2 or _ => new I3(0, hsize, 0),
        };
        var edgeDirs = dir switch {
            0 => (1, 2),
            1 => (2, 0),
            2 or _ => (0, 1),
        };
        outArray[0] = new EdgeCoord(center - axis1, hsize, edgeDirs.Item1);
        outArray[1] = new EdgeCoord(center, hsize, edgeDirs.Item1);
        outArray[2] = new EdgeCoord(center - axis2, hsize, edgeDirs.Item2);
        outArray[3] = new EdgeCoord(center, hsize, edgeDirs.Item2);
    }

    // subdidivide a face into 4 children
    private void SubdivideFace(int faceUid, I3 cellMin, int cellSize, int faceIndex) {
        Debug.Assert(cellSize > 1 && (cellSize&1) == 0);
        // NOTE: assuming all face children are allocated continuously
        int firstFaceChild = newFace();
        int c2 = newFace();
        Debug.Assert(c2 == firstFaceChild + 1);
        int c3 = newFace();
        Debug.Assert(c3 == firstFaceChild + 2);
        int c4 = newFace();
        Debug.Assert(c4 == firstFaceChild + 3);
        var faceSpan = CollectionsMarshal.AsSpan(facePool);
        var faceEdgeSpan = CollectionsMarshal.AsSpan(faceEdgePool);
        ref FaceEdges face = ref faceEdgeSpan[faceUid];
        ref FaceEdges face0 = ref faceEdgeSpan[firstFaceChild + 0];
        ref FaceEdges face1 = ref faceEdgeSpan[firstFaceChild + 1];
        ref FaceEdges face2 = ref faceEdgeSpan[firstFaceChild + 2];
        ref FaceEdges face3 = ref faceEdgeSpan[firstFaceChild + 3];

        var hsize = cellSize >> 1;

        var newEdges = new EdgeCoord[4];
        FaceInternalEdges(faceIndex, cellSize, cellMin, newEdges);

        // inner
        face0[2] = newEdges[0];
        face0[3] = newEdges[2];
        face1[1] = newEdges[2];
        face1[2] = newEdges[1];
        face2[0] = newEdges[0];
        face2[3] = newEdges[3];
        face3[1] = newEdges[3];
        face3[0] = newEdges[1];

        // outer
        var (e01, e02) = face[0].Subdivide();
        face0[0] = e01;
        face1[0] = e02;
        var (e11, e12) = face[1].Subdivide();
        face0[1] = e11;
        face2[1] = e12;
        var (e21, e22) = face[2].Subdivide();
        face2[2] = e21;
        face3[2] = e22;
        var (e31, e32) = face[3].Subdivide();
        face1[3] = e31;
        face3[3] = e32;

        faceSpan[faceUid].firstChild = firstFaceChild;
    }

    // collect all segments from the quadtree
    private void GetCellFaceSegments(int cellId, List<Segment> outSegments, List<Segment>[] faceSegments) {
        bool reverse = false;
		void recurseFace(int faceId) {
			if (facePool[faceId].IsLeaf) {
                var segs = faceSegments[faceId];
                var len = segs.Count;
                // segments on negative face need to be reversed
                for (int i = 0; i < len; i++) {
                    var cseg = reverse ? segs[i] : segs[len - i - 1];
                    if (reverse) {
                        outSegments.Add(cseg);
                    } else {
                        // reverse segment
                        outSegments.Add(new Segment(cseg.v2, cseg.v1, cseg.nDir));
                    }
                }
			} else {
                int fc = facePool[faceId].firstChild;
                recurseFace(fc);
                recurseFace(fc + 1);
                recurseFace(fc + 2);
                recurseFace(fc + 3);
			}
		}
        for (int i = 0; i < 6; i++) {
            var face = cellFacePool[cellId][i];
            reverse = (i&1) == 0;
			recurseFace(face);
		}
    }

    // finds all segment loops in CMS cell
    private static Span<int[]> GetLoops(in CellCtx ctx) {
        var ret = new int[4][];
        int curRetIdx = 0;
        var visited = ctx.visited;
        var segments = ctx.cellSegments;
        visited.Clear();
        ctx.vertices.Clear();
		foreach(var seg in segments) {
            ctx.vertices.Add(seg.v1);
            ctx.vertices.Add(seg.v2);
		}

		static int[] BuildLoop(in CellCtx ctx, Vertex start)
		{
			var loop = ctx.loop;
            loop.Clear();
            var usedSegments = ctx.usedSegments;
			Vertex current = start;
			Vertex? previous = null;

			while (true)
			{
				ctx.visited.Add(current);

                Segment? nextSeg = null;
                foreach(var seg in ctx.cellSegments) {
                    if ((seg.v1.Equals(current) || seg.v2.Equals(current))
                        && !usedSegments.Contains(seg))
                    {
                        nextSeg = seg;
                        break;
                    }
                }

				if (nextSeg == null) {
                    //throw new Exception("No next segment in loop");
					break; // dead end or loop complete
				}

				usedSegments.Add(nextSeg.Value);
				loop.Add(ctx.cellSegments.IndexOf(nextSeg.Value));

				Vertex next = nextSeg.Value.v1.Equals(current) ? nextSeg.Value.v2 : nextSeg.Value.v1;

				if (next.Equals(start))
				{
					// closed loop
					break;
				}

				previous = current;
				current = next;
			}
			return loop.ToArray();
		}
        ctx.usedSegments.Clear();

		foreach (var v in ctx.vertices)
		{
			if (!visited.Contains(v))
			{
				var loop = BuildLoop(ctx, v);
                ret[curRetIdx++] = loop;
			}
		}
#if DEBUG
        Debug.Assert(ctx.usedSegments.Count == segments.Count, "Not all segments were used");
#endif
		return ret.AsSpan().Slice(0, curRetIdx);
	}
}
