namespace Aviz.Cms;

using V3 = System.Numerics.Vector3;
using V2 = System.Numerics.Vector2;

internal static class Util
{
    // Möller–Trumbore intersection algorithm
    public static bool SegmentTriangleIntersect(V3 p0, V3 p1, V3 v0, V3 v1, V3 v2)
    {
        const float eps = 1e-6f;

        V3 dir =  p1 - p0;
        V3 edge1 =  v1 - v0;
        V3 edge2 =  v2 - v0;

        V3 h = V3.Cross(dir, edge2);
        float a = V3.Dot(edge1, h);

        if (a > -eps && a < eps) return false; // parallel

        float f = 1.0f / a;
        V3 s =  p0 - v0;
        float u = f * V3.Dot(s, h);
        if (u < 0.0f || u > 1.0f) return false;

        V3 q = V3.Cross(s, edge1);
        float v = f * V3.Dot(dir, q);
        if (v < 0.0f || u + v > 1.0f) return false;

        float t = f * V3.Dot(edge2, q);
        if (t < 0.0f || t > 1.0f) return false; // outside segment

        return true;
    }

    // NOTE: this doesn't handle coplanar triangles
    //   (which seems fine for our use case)
    public static bool TriTriIntersect(V3 v0, V3 v1, V3 v2, V3 u0, V3 u1, V3 u2)
    {
        // edges of first triangle
        if (SegmentTriangleIntersect(v0, v1, u0, u1, u2)) return true;
        if (SegmentTriangleIntersect(v1, v2, u0, u1, u2)) return true;
        if (SegmentTriangleIntersect(v2, v0, u0, u1, u2)) return true;
        // Edges of second triangle
        if (SegmentTriangleIntersect(u0, u1, v0, v1, v2)) return true;
        if (SegmentTriangleIntersect(u1, u2, v0, v1, v2)) return true;
        if (SegmentTriangleIntersect(u2, u0, v0, v1, v2)) return true;
        return false;
    }

    // Get run lengths, for example
    //   333322111 -> 4,2,3
    public static IEnumerable<int> RunLengths<T>(IEnumerable<T> seq) {
        using var e = seq.GetEnumerator();
        if (!e.MoveNext())
            yield break;

        var comparer = EqualityComparer<T>.Default;
        T prev = e.Current;
        int count = 1;

        while (e.MoveNext()) {
            if (comparer.Equals(e.Current, prev)) {
                count++;
            } else {
                yield return count;
                count = 1;
                prev = e.Current;
            }
        }
        yield return count;
    }

    public static bool SegmentsIntersect(V2 s1, V2 e1, V2 s2, V2 e2)
    {
        var d1 = e1 - s1;
        var d2 = e2 - s2;
        var d = s1 - s2;
        var denom = d2.Y * d1.X - d2.X * d1.Y;
        if (denom == 0.0f)
        {
            return false;
        }
        var t1 = (d2.X * d.Y - d2.Y * d.X) / denom;
        var t2 = (d1.X * d.Y - d1.Y * d.X) / denom;
        return t1 >= 0.0f && t1 <= 1.0f && t2 >= 0.0f && t2 <= 1.0f;
    }
}
