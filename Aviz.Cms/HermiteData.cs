using System;

namespace Aviz.Cms;

using V3 = System.Numerics.Vector3;
using Debug = System.Diagnostics.Debug;

public class HermiteData {
    // NOTE: might be better to use Dictionary for very sparse geometry
	public int[,,,] intersections;
    //Dictionary<(int, int, int, int), (V3 p, V3 n)> intersections = new();
    List<(V3 p, V3 n)> intersectionList = new();
	public bool[,,] isInside;
	public V3 offset;
	public float step;
    public int size;

    public HermiteData(bool[,,] isInside, V3 offset, float step) {
        this.isInside = isInside;
        this.offset = offset;
        this.step = step;
        this.size = isInside.GetLength(0)-1;
        this.intersections = new int[3, size+1, size+1, size];
        // dummy intersection at index 0
        // avoids setting default value to -1
        this.intersectionList.Add((V3.Zero, V3.Zero));
    }

	public static HermiteData FromSdf(Func<V3, float> eval, Func<V3, V3> evalNormal, V3 offset, float step, int gridSize) {
		var gp1 = gridSize+1;
        var isInside = new bool[gp1, gp1, gp1];
        var ret = new HermiteData(isInside, offset, step);

		// eval grid
		for (int i = 0; i < gp1; i++)
		for (int j = 0; j < gp1; j++)
		for (int k = 0; k < gp1; k++) {
			var location = new V3(offset.X + i * step, offset.Y + j * step, offset.Z + k * step);
			ret.isInside[i, j, k] = eval(location) <= 0.0f;
		}

		for (int dir = 0; dir < 3; dir++) {
			for (int k = 0; k < gridSize; k++)
			for (int i = 0; i < gp1; i++) {
				for (int j = 0; j < gp1; j++) {
					var sample1 = dir switch {
						0 => ret.isInside[k, i, j],
						1 => ret.isInside[i, k, j],
						2 => ret.isInside[i, j, k],
						_ => false
					};
					var sample2 = dir switch {
						0 => ret.isInside[k+1, i, j],
						1 => ret.isInside[i, k+1, j],
						2 => ret.isInside[i, j, k+1],
						_ => false
					};
					if (sample1 == sample2) continue;
					var start = dir switch {
						0 => new V3(offset.X + k * step, offset.Y + i * step, offset.Z + j * step),
						1 => new V3(offset.X + i * step, offset.Y + k * step, offset.Z + j * step),
						2 or _ => new V3(offset.X + i * step, offset.Y + j * step, offset.Z + k * step),
					};
					var end = start + (dir switch {
						0 => new V3(step, 0.0f, 0.0f),
						1 => new V3(0.0f, step, 0.0f),
						2 or _ => new V3(0.0f, 0.0f, step),
					});
					float curMin = float.MaxValue;
					V3 minPos = V3.Zero;
					for (int s = 0; s <= 100; s++) {
						V3 pos = V3.Lerp(start, end, s / 100.0f);
						float val = eval(pos);
						float dist = val*val;
						if (dist < curMin) {
							curMin = dist;
							minPos = pos;
						}
					}
					V3 minNormal = evalNormal(minPos);
                    var newIndex = ret.intersectionList.Count;
                    ret.intersectionList.Add((minPos, minNormal));
					ret.intersections[dir, i, j, k] = newIndex;
				}
			}
		}
        return ret;
	}

    public (V3 p, V3 n)? GetIntersection(EdgeCoord coord) {
        for (int i = 0; i < coord.count; i++) {
            int curI = coord.dir switch {
                0 => this.intersections[0, coord.y, coord.z, coord.x+i],
                1 => this.intersections[1, coord.x, coord.z, coord.y+i],
                2 or _ => this.intersections[2, coord.x, coord.y, coord.z+i],
            };
            if (curI > 0) {
                return this.intersectionList[curI];
            }
        }
        return null;
    } 

    public (EdgeCoord leaf, float t, V3 n)? GetIntersectionWithLeaf(EdgeCoord coord) {
        for (int i = 0; i < coord.count; i++) {
            int curI = coord.dir switch {
                0 => this.intersections[0, coord.y, coord.z, coord.x+i],
                1 => this.intersections[1, coord.x, coord.z, coord.y+i],
                2 or _ => this.intersections[2, coord.x, coord.y, coord.z+i],
            };
            var leaf = coord.dir switch {
                0 => new EdgeCoord() { x = coord.x+i, y = coord.y, z = coord.z, count = 1, dir = coord.dir },
                1 => new EdgeCoord() { x = coord.x, y = coord.y+i, z = coord.z, count = 1, dir = coord.dir },
                2 or _ => new EdgeCoord() { x = coord.x, y = coord.y, z = coord.z+i, count = 1, dir = coord.dir },
            };
            if (curI > 0) {
                var cur = this.intersectionList[curI];
                V3 start  = offset + new V3(leaf.x*step, leaf.y*step, leaf.z*step);
                float t = (cur.p - start).Length() / step;
                return (leaf, t, cur.n);
            }
        }
        return null;
    }

    // get intersection in cell-local coordinates
    public (V3 p, V3 n)? GetIntersectionNormalized(EdgeCoord coord, I3 cellMin, int cellSize) {
        for (int i = 0; i < coord.count; i++) {
            int curI = coord.dir switch {
                0 => this.intersections[0, coord.y, coord.z, coord.x+i],
                1 => this.intersections[1, coord.x, coord.z, coord.y+i],
                2 or _ => this.intersections[2, coord.x, coord.y, coord.z+i],
            };
            if (curI > 0) {
                var curv = this.intersectionList[curI];
                V3 normPos = curv.p - offset - cellMin.ToV3()*step;
                Debug.Assert(float.Abs(cellSize) > 1e-6);
                normPos = normPos / ((float)cellSize * step);
                return (normPos, curv.n);
            }
        }
        return null;
    }

    // convert normalized vertex back to global position
    public V3 ResolveVertex(NewVertex vert) {
        return vert.cellMin.ToV3()*step + offset + vert.position*vert.cellSize*step;
    }

    internal (float t, V3 n)[] FindEdgeIntersections(EdgeCoord coord) {
        Debug.Assert(coord.count > 0);
        List<(float t, V3 n)> intersectionsList = new();
        for (int i = 0; i < coord.count; i++) {
            int curI = coord.dir switch {
                0 => this.intersections[0, coord.y, coord.z, coord.x + i],
                1 => this.intersections[1, coord.x, coord.z, coord.y + i],
                2 or _ => this.intersections[2, coord.x, coord.y, coord.z + i],
            };
            if (curI > 0) {
                var cur = this.intersectionList[curI];
                V3 start  = offset + new V3(coord.x*step, coord.y*step, coord.z*step);
                float t = (cur.p - start).Length() / (coord.count*step);
                intersectionsList.Add((t, cur.n));
            }
        }
        return intersectionsList.ToArray();
    }
}
