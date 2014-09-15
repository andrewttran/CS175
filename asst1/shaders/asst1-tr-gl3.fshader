#version 130

//uniform float uVertexScale;
uniform sampler2D uTex2;

in vec2 vTexCoord;
in vec2 vTemp;
//in vec3 aColor;

out vec4 fragColor;

void main(void) {

  vec4 texColor2 = texture(uTex2, vTexCoord);
  vec4 someColor = vec4(vTexCoord.x, vTexCoord.y, 0.5, 1);
  fragColor = texColor2 + someColor;

}
