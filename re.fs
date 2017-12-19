#version 120


uniform sampler2D QueryFluxRadiusTexture;
uniform sampler2D QueryEmissionPhotonCountTexture;
uniform float TotalPhotonNum;


float sRGB(const float c)
{
	if (c < 0.0031308)
	{
		return 12.92 * c;
	}
	else
	{
		const float a = 0.055;
		return (1.0 + a) * pow(c, 1.0 / 2.4) - a;
	}
}


void main()
{
	vec2 PixelIndex = gl_TexCoord[0].st;

	// fetch various data of the measurement point
	vec4 QueryFluxRadius = texture2D(QueryFluxRadiusTexture, PixelIndex);
	float QueryRadius = QueryFluxRadius.w;
	vec3 QueryFlux = QueryFluxRadius.xyz;
	vec3 QueryEmission = texture2D(QueryEmissionPhotonCountTexture, PixelIndex).rgb;

	// perform progressive density estimation
	gl_FragColor = vec4(QueryFlux / (QueryRadius * QueryRadius * 3.141592 * TotalPhotonNum), 1.0);

	// add emission
	gl_FragColor = gl_FragColor + vec4(QueryEmission, 0.0); 

	// tone mapping
	const float Exposure = 60000.0;
	gl_FragColor = vec4(1.0) - exp(-gl_FragColor * Exposure);

	// sRGB conversion
	gl_FragColor.r = sRGB(gl_FragColor.r);
	gl_FragColor.g = sRGB(gl_FragColor.g);
	gl_FragColor.b = sRGB(gl_FragColor.b);
}
