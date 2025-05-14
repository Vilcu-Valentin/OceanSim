using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter))]
public class ProceduralOcean : MonoBehaviour
{
	private Mesh mesh;
	private Vector3[] baseVertices;
	private Vector3[] displacedVertices;

	public float waveHeight = 1f;
	public float waveFrequency = 1f;
	public float waveSpeed = 1f;

	void Start()
	{
		mesh = GetComponent<MeshFilter>().mesh;
		baseVertices = mesh.vertices;
		displacedVertices = new Vector3[baseVertices.Length];
	}

	void Update()
	{
		for (int i = 0; i < baseVertices.Length; i++)
		{
			Vector3 vertex = baseVertices[i];
			vertex.y = Mathf.Sin(Time.time * waveSpeed + vertex.x * waveFrequency + vertex.z * waveFrequency) * waveHeight;
			displacedVertices[i] = vertex;
		}

		mesh.vertices = displacedVertices;
		mesh.RecalculateNormals(); // pentru iluminare corectă

		GetComponent<MeshCollider>().sharedMesh = null;
		GetComponent<MeshCollider>().sharedMesh = mesh;
	}
}