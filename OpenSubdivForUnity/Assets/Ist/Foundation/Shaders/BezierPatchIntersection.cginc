#ifndef BezierPatchIntersection_h
#define BezierPatchIntersection_h

// thanks to @ototoi, @gishicho
// http://jcgt.org/published/0004/01/04/

// prototypes

struct BezierPatchHit {
    float t, u, v;
    int clip_level;
};

bool BPIRaycast(BezierPatch bp, Ray ray, real zmin, real zmax, real eps, out BezierPatchHit hit);




// implements

#ifndef BPI_MAX_STACK_DEPTH
    #define BPI_MAX_STACK_DEPTH 40
#endif
#ifndef BPI_MAX_LOOP
    #define BPI_MAX_LOOP 2000
#endif

struct BPIWorkingBuffer {
    BezierPatch source; // input
    BezierPatch crop;
    BezierPatch rotate;
    float4 uv_range; // input
};


float3x3 BPIRotate2D_(float3 dx) {
    float2 x = normalize(dx.xy);
    float2 y = float2(-x[1], x[0]);
    return float3x3(
        x[0], y[0], 0.0,
        x[1], y[1], 0.0,
         0.0,  0.0, 1.0
    );
}

#if defined(_INTERSECTION_BILINEAR)

inline int BPISolve2_(out real2 root, real3 coeff) {
	real A = coeff.x;
	real B = coeff.y;
	real C = coeff.z;
	if (abs(A) <= 1e-6) {
		if (abs(B) <= 1e-6) return 0;
		real x = -C / B;
		root = real2(x, 0);
		return 1;
	} else {
		real D = B * B - 4 * A * C;
		if (D < 0) {
			return 0;
		} else if (D == 0) {
			real x = -0.5 * B / A;
			root = real2(x, 0);
			return 1;
		} else {
			real x1 = (abs(B) + sqrt(D)) / (2.0 * A);
			if (B >= 0) {
				x1 = -x1;
			}
			real x2 = C / (A * x1);
			root = real2(min(x1, x2), max(x1, x2));
			return 2;
		}
	}
}

inline real BPIComputeU_(real A1, real A2, real B1, real B2, real C1, real C2, real D1, real D2, real v) {
	//return div((v*(C1-C2)+(D1-D2)),(v*(A2-A1)+(B2-B1)));
	real a = v * A2 + B2;
	real b = v * (A2 - A1) + B2 - B1;
	return (abs(b) >= abs(a))
		? (v * (C1 - C2) + D1 - D2) / b
		: (-v * C2 - D2) / a;
}

inline real BPIComputeT_(real a, real b, real c, real d, real iq, real2 uv) {
	return ((uv.x * uv.y) * a + uv.x * b + uv.y * c + d) * iq;
}

inline bool BPISolveBilinearPatch_(out real3 uvt, real tmin, real tmax, real3 uuvvtt) {
	real eps = 1e-5;
	if (tmin - eps <= uuvvtt.z && uuvvtt.z <= tmax + eps) {
		uvt = uuvvtt;
		return true;
	}
	uvt = real3(0, 0, 0);
	return false;
}

inline bool BPITestBilinearPatch_(out real3 uvt, real3 p[4], real tmin, real tmax, real uvMargin) {
	uvt = real3(0, 0, 0);
	const real3 p00 = p[0];
	const real3 p10 = p[1];
	const real3 p01 = p[2];
	const real3 p11 = p[3];

	static const int nPlane = 2;
	real3 a = p11 - p10 - p01 + p00;
	real3 b = p10 - p00;
	real3 c = p01 - p00;
	real3 d = p00;

	//xz-zx
	real A1 = a[0];
	real B1 = b[0];
	real C1 = c[0];
	real D1 = d[0];

	//yz-zy
	real A2 = a[1];
	real B2 = b[1];
	real C2 = c[1];
	real D2 = d[1];

	real F1 = A2 * C1 - A1 * C2;
	real F2 = A2 * D1 - A1 * D2 + B2 * C1 - B1 * C2;
	real F3 = B2 * D1 - B1 * D2;

	real2 root;
	int count = BPISolve2_(root, real3(F1, F2, F3));

	if (count <= 0) return false;

	bool hit = false;
	real3 uuvvtt = real3(0, 0, 0);
	[unroll(2)]
	for (int i = 0; i < count; ++i) {
		uuvvtt.y = root[i];
		if (0 - uvMargin <= uuvvtt.y && uuvvtt.y <= 1 + uvMargin) {//TODO
			uuvvtt.y = saturate(uuvvtt.y);
			uuvvtt.x = BPIComputeU_(A1, A2, B1, B2, C1, C2, D1, D2, uuvvtt.y);
			if (0 - uvMargin <= uuvvtt.x && uuvvtt.x <= 1 + uvMargin) {//TODO
				uuvvtt.x = saturate(uuvvtt.x);
				uuvvtt.z = BPIComputeT_(a[nPlane], b[nPlane], c[nPlane], d[nPlane],
					1.0, uuvvtt.xy);
				if (BPISolveBilinearPatch_(uvt, tmin, tmax, uuvvtt)) {
					tmax = uvt.z;
					hit = true;
				}
			}
		}
	}
	return hit;
}

inline bool BPITestBezierClipL_(
	BezierPatch bp, out real3 uvt, real2 uv0, real2 uv1, real zmin, real zmax) {
	real3 p[4];
	real3 ray_org = real3(0.0, 0.0, 0.0);
	real3 ray_dir = real3(0.0, 0.0, 1.0);
	p[0] = bp.cp[0];
	p[1] = bp.cp[3];
	p[2] = bp.cp[12];
	p[3] = bp.cp[15];
	real m_uvMargin = 1e-1; // ???
	if (BPITestBilinearPatch_(uvt, p, zmin, zmax, m_uvMargin)) {
		uvt.xy = lerp(uv0, uv1, uvt.xy);
		return true;
	}
	uvt = real3(0, 0, zmax);
	return false;
}

#else

inline bool BPITriangleIntersect_(
	inout real tout, out real uout, out real vout,
	real3 p0, real3 p1, real3 p2,
	real3 ray_org, real3 ray_dir) {
	real3 e1, e2;
	real3 p, s, q;

	e1 = p1 - p0;
	e2 = p2 - p0;
	p = cross(ray_dir, e2);

	real det = dot(e1, p);
	real inv_det = 1.0 / det;

	s = ray_org - p0;
	q = cross(s, e1);

	real u = dot(s, p) * inv_det;
	real v = dot(q, ray_dir) * inv_det;
	real t = dot(e2, q) * inv_det;

	real eps = 1e-2;
	if (u < 0.0 - eps || u > 1.0 + eps) return false;
	if (v < 0.0 - eps || u + v > 1.0 + eps) return false;
	if (t < 0.0 - eps || t > tout + eps) return false;

	tout = t;
	uout = u;
	vout = v;
	return true;
}

inline bool BPITestBezierClipL_(
    BezierPatch bp, out real3 uvt, real2 uv0, real2 uv1, real zmin, real zmax) {
	real3 p0, p1, p2, p3;
	real3 ray_org = real3(0.0, 0.0, 0.0);
	real3 ray_dir = real3(0.0, 0.0, 1.0);
    p0 = bp.cp[0];
    p1 = bp.cp[3];
    p2 = bp.cp[12];
    p3 = bp.cp[15];
    bool ret = false;
	real t = zmax, uu = 0.0, vv = 0.0;
    if (BPITriangleIntersect_(t, uu, vv, p0, p2, p1, ray_org, ray_dir)) {
		real ww = 1.0 - (uu + vv);
		real u = ww*0.0 + uu*0.0 + vv*1.0; //00 - 01 - 10
		real v = ww*0.0 + uu*1.0 + vv*0.0; //00 - 01 - 10
		uvt.xy = lerp(uv0.xy, uv1.xy, real2(u, v));
		uvt.z = t;
        ret = true;
    }
    if (BPITriangleIntersect_(t, uu, vv, p1, p2, p3, ray_org, ray_dir)) {
		real ww = 1.0 - (uu + vv);
		real u = ww*1.0 + uu*0.0 + vv*1.0; //10 - 01 - 11
		real v = ww*0.0 + uu*1.0 + vv*1.0; //10 - 01 - 11
        uvt.xy = lerp(uv0.xy, uv1.xy, real2(u, v));
        uvt.z = t;
        ret = true;
    }
    return ret;
}
#endif

inline bool BPITestBounds_(inout BPIWorkingBuffer work, inout BezierPatchHit info, real zmin, real zmax, real eps) {
	real3 bmin, bmax;
    BPGetMinMax(work.source, bmin, bmax, eps * 1e-3);
	return !(0.0 < bmin.x || bmax.x < 0.0 || 0.0 < bmin.y || bmax.y < 0.0 || bmax.z < zmin || zmax < bmin.z);
}

inline bool BPITestBezierPatch_(inout BPIWorkingBuffer work, inout BezierPatchHit info, real zmin, real zmax, real eps) {
    info = (BezierPatchHit)0;

    // non-recursive iteration
    float4 range_stack[BPI_MAX_STACK_DEPTH];
    int stack_index = 0;
    range_stack[0] = work.uv_range;

    bool ret = false;
	[loop]
    for (int i = 0; i < BPI_MAX_LOOP && stack_index >= 0; ++i) {

        // pop a patch range and crop
        float u0 = range_stack[stack_index].x;
        float u1 = range_stack[stack_index].y;
        float v0 = range_stack[stack_index].z;
        float v1 = range_stack[stack_index].w;
        --stack_index;

        BPCrop(work.source, work.crop, float2(u0, v0), float2(u1, v1));
        float3 LU = work.crop.cp[3] - work.crop.cp[0];
        float3 LV = work.crop.cp[12] - work.crop.cp[0];
        bool clipU = length(LU) > length(LV);

        float3 bmin, bmax;
        // rotate and bmin/bmax
        float3 dx = clipU
            ? work.crop.cp[12] - work.crop.cp[0] + work.crop.cp[15] - work.crop.cp[3]
            : work.crop.cp[3] - work.crop.cp[0] + work.crop.cp[15] - work.crop.cp[12];
        work.rotate = work.crop;
        BPTransform(work.rotate, BPIRotate2D_(dx));
        BPGetMinMax(work.rotate, bmin, bmax, eps*1e-3);

        // out
        if (0.0 < bmin.x || bmax.x < 0.0 || 0.0 < bmin.y || bmax.y < 0.0 || bmax.z < zmin || zmax < bmin.z) {
            continue;
        }

        // if it's small enough, test bilinear.
        if ((bmax.x - bmin.x) < eps || (bmax.y - bmin.y) < eps) {
			real3 uvt;
            if (BPITestBezierClipL_(work.crop, uvt, real2(u0, v0), real2(u1, v1), zmin, zmax)) {
                info.u = uvt.x;
                info.v = uvt.y;
                info.t = uvt.z;
                zmax = info.t;
                ret = true;
            }
            // find another intersection
            continue;
        }
        info.clip_level = i;

        // push children ranges
        if (clipU) {
            float um = (u0 + u1) * 0.5;
            range_stack[++stack_index] = float4(u0, um, v0, v1);
            range_stack[++stack_index] = float4(um, u1, v0, v1);
        } else {
            float vm = (v0 + v1) * 0.5;
            range_stack[++stack_index] = float4(u0, u1, v0, vm);
            range_stack[++stack_index] = float4(u0, u1, vm, v1);
        }

        if (stack_index >= BPI_MAX_STACK_DEPTH - 1) break;
    }
    return ret;
}

inline bool BPIRaycast(BezierPatch bp, Ray ray, real zmin, real zmax, real eps, out BezierPatchHit hit) {
    BPIWorkingBuffer work = (BPIWorkingBuffer)0;

    work.source = bp;
    work.uv_range = float4(0.0, 1.0, 0.0, 1.0);
    BPTransform(work.source, ZAlign(ray.origin, ray.direction));

    //// all pixels pass this test when draw aabb as mesh
    //// (but viable if run on GLSLSandbox etc.)
    //if (!BPITestBounds_(work, hit, zmin, zmax, BPI_EPS)) {
    //    return false;
    //}

    hit.t = zmax;
	return BPITestBezierPatch_(work, hit, zmin, zmax, eps);
}

#endif // BezierPatchIntersection_h
