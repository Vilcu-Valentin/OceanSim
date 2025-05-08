using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FloatingObject : MonoBehaviour
{
    public float offset = 0.5f;
    public Transform[] floatPoints; // Cele 4 puncte pentru Raycast
    public LayerMask waterLayer; // Layer-ul la care este oceanul

    void Update()
    {
        if (floatPoints == null || floatPoints.Length == 0)
            return;

        Vector3 averageHitPoint = Vector3.zero;
		Vector3[] hitPoints = new Vector3[4];
		int hitCount = 0;

        foreach(Transform point in floatPoints)
        {
            if (Physics.Raycast(point.position, Vector3.down, out RaycastHit hit, 100f, waterLayer))
            {
                averageHitPoint += hit.point;
                hitCount++;

				Debug.Log("Hit water at: " + hit.point);
				Debug.DrawRay(point.position, Vector3.down * hit.distance, Color.red);
            }
        }

        if (hitCount > 0)
        {
            averageHitPoint /= hitCount;

            // Pozitioneaza barca la media pozitiilor atinse
            Vector3 newPos = new Vector3(transform.position.x, averageHitPoint.y + offset, transform.position.z);
            transform.position = Vector3.Lerp(transform.position, newPos, Time.deltaTime * 5f);

			// Rotatie
			Vector3 forward = (hitPoints[2] + hitPoints[3]) * 0.5f - (hitPoints[0] + hitPoints[1]) * 0.5f;
			Vector3 right = (hitPoints[1] + hitPoints[3]) * 0.5f - (hitPoints[0] + hitPoints[2]) * 0.5f;

			Vector3 up = Vector3.Cross(right, forward).normalized;

			Quaternion targetRotation = Quaternion.LookRotation(forward, up);
			transform.rotation = Quaternion.Slerp(transform.rotation, targetRotation, Time.deltaTime * 2f);
		}
    }
}
