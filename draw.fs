#version 120


uniform sampler2D input_tex;


void main()
{
	gl_FragColor = texture2D(input_tex, gl_TexCoord[0].st) / 1.0;
}
