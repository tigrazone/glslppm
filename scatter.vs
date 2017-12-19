#version 120


uniform sampler2D PhotonIndexTexture;
uniform vec4 BufInfo;

uniform float PhotonBufferSize;
varying vec4 p;


void main()
{
	// read hashed photon index
	vec2 TexCoord = (gl_Vertex.xy + vec2(1.0)) * 0.5;
	vec4 PhotonIndex = texture2D(PhotonIndexTexture, TexCoord);
	vec2 PhotonListIndex = PhotonIndex.zw; 

	// global 1d index in the photon buffer (i.e., photon id)
	float z = (PhotonIndex.x / PhotonBufferSize) + PhotonIndex.y;

	gl_Position = vec4(PhotonListIndex * BufInfo.zw * 2.0 - vec2(1.0), z, 1.0);
	p = PhotonIndex;
}
