#version 130

//uniform float uVertexScale;
uniform sampler2D uTex2;

in vec2 vTexCoord;
in vec2 vTemp;

out vec4 fragColor;

void main(void) {
  vec4 color = vec4(vTemp.x, vTemp.y, 0.5, 1);

  // fragColor is a vec4. The components are interpreted as red, green, blue, and alpha
  fragColor = texture(uTex2, vTexCoord);
}
