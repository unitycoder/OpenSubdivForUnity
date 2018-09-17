using UnityEngine;

namespace Ist {
	[ExecuteInEditMode]
	public class RandomBezierPatch : IBezierPatchContainer {
		[SerializeField] BezierPatch[] patches = new BezierPatch[0];

		BezierPatchRaw[] m_bpatch_raw;
		BezierPatchAABB[] m_aabb;
		bool m_dirty = true;

		[Range(1, 256)]
		public int Amount = 16;

		Vector3[,] velocities;

		private void OnValidate() {
			m_bpatch_raw = null;
		}

		void GenerateAtRandom() {
			patches = new BezierPatch[Amount];
			velocities = new Vector3[Amount, 16];
			for (int i = 0; i < Amount; i++) {
				if (patches[i] == null)
					patches[i] = new BezierPatch();
				for (int j = 0; j < 16; j++) {
					patches[i].cp[j] = Random.insideUnitSphere;
					velocities[i, j] = Random.onUnitSphere;
				}
			}
		}

		void UpdateBezierPatch() {
			if (!m_dirty) { return; }
			if (m_bpatch_raw == null || m_bpatch_raw.Length != Amount) {
				m_bpatch_raw = new BezierPatchRaw[Amount];
				m_aabb = new BezierPatchAABB[Amount];
				GenerateAtRandom();
			}

			var localToWorldMatrix = GetComponent<Transform>().localToWorldMatrix;
			for (int i = 0; i < Amount; i++) {
				patches[i].GetRawData(ref m_bpatch_raw[i]);
				m_bpatch_raw[i].Transform(ref localToWorldMatrix);
				m_bpatch_raw[i].GetAABB(ref m_aabb[i]);
			}
			m_dirty = false;
		}

		public override BezierPatchRaw[] GetBezierPatches() {
			UpdateBezierPatch();
			return m_bpatch_raw;
		}

		public override BezierPatchAABB[] GetAABBs() {
			UpdateBezierPatch();
			return m_aabb;
		}


		public BezierPatch[] Patches { get { return patches; } }

		void Update() {
			m_dirty = true;
			UpdateBezierPatch();
			var t = Time.deltaTime;
			for (int i = 0; i < Amount; i++) {
				for (int j = 0; j < 16; j++) {
					var p = patches[i].cp[j];
					if (p.sqrMagnitude >= 1 && (p + velocities[i, j]).sqrMagnitude >= 1) {
						velocities[i, j] = -velocities[i, j];
					}

					p += t * velocities[i, j];
					patches[i].cp[j] = p;
				}
			}
		}

	}
}
