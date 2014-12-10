#version 130

uniform float uWidthScale, uHeightScale, uShiftX, uShiftY;

in vec2 aPosition;
in vec2 aTexCoord;

out vec2 vTexCoord;
out vec2 vTemp;

void main() {
  gl_Position = vec4(aPosition.x * uWidthScale + uShiftX, aPosition.y * uHeightScale + uShiftY, 0, 1);
  vTexCoord = aTexCoord;
  vTemp = vec2(1, 1);
}
