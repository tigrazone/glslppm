#version 120


uniform sampler2D PhotonIndexTexture;
uniform vec4 BufInfo;


void main()
{
	vec2 TexCoord = (gl_Vertex.xy + vec2(1.0)) * 0.5;
	vec2 PhotonListIndex = texture2D(PhotonIndexTexture, TexCoord).zw;
	gl_Position = vec4(PhotonListIndex * BufInfo.zw * 2.0 - vec2(1.0), 0.5, 1.0);
}
