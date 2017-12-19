#version 120


uniform sampler2D Texture;
uniform vec2 Offset;


void main()
{
	vec2 t0 = gl_TexCoord[0].st; 
	vec2 t1 = gl_TexCoord[0].st + vec2(Offset.x, 0.0);
	vec2 t2 = gl_TexCoord[0].st + vec2(0.0, Offset.y);
	vec2 t3 = gl_TexCoord[0].st + Offset;

	vec4 v0 = texture2D(Texture, t0);
	vec4 v1 = texture2D(Texture, t1);
	vec4 v2 = texture2D(Texture, t2);
	vec4 v3 = texture2D(Texture, t3);

	gl_FragColor = max(max(v0, v1), max(v2, v3));
}
