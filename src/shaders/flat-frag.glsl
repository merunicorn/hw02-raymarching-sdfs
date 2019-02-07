#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

vec3 rayCast(vec4 s) {
    // multiply by far clip plane
    float far_clip = 1000.0;
    float near_clip = 0.1;
    s *= far_clip;

    // multiply by inverse projection matrix
    mat4 proj;
    float fov = 90.0;
    float fov_rad = (fov * 3.14159) / 180.0;
    float S_f = tan(fov_rad / 2.0);
    float a = 1.0 / ((u_Dimensions.x / u_Dimensions.y) * S_f);
    float b = 1.0 / S_f;
    float P = far_clip / (far_clip - near_clip);
    float Q = (-1.f * far_clip * near_clip) / (far_clip - near_clip);
    proj[0][0] = a;
    proj[1][1] = b;
    proj[2][2] = P;
    proj[3][2] = Q;
    proj[2][3] = 1.0;
    proj[3][3] = 0.0;
    s = inverse(proj) * s;

    // multiply by inverse of view matrix
    mat4 view;
    mat4 orient;
    mat4 transl;
    vec3 forw_axis = u_Ref - u_Eye;
    vec3 right_axis = cross(u_Up, forw_axis);
    vec3 up_axis = cross(right_axis, forw_axis);

    // special case where forward is world up
    if (forw_axis == u_Up) {
      right_axis = vec3(1.0, 0.0, 0.0);
      up_axis = vec3(0.0, 0.0, -1.0);
    }

    // set up orientation and translation matrices
    for (int i = 0; i < 3; i++) {
      orient[0][i] = right_axis[i];
      orient[1][i] = up_axis[i];
      orient[2][i] = forw_axis[i];
      transl[3][i] = u_Eye[i] * -1.0;
    }
    view = orient * transl;
    s = inverse(view) * s;

    //set up ray
    vec3 origin = u_Eye;
    vec3 dir = normalize(vec3(s) - u_Eye);
    return dir;
}

float sdf_sphere(vec3 p, float rad) {
  float sphere = length(p) - rad;
  return sphere;
}

float sdf_box(vec3 p, vec3 b) {
  vec3 d = abs(p) - b;
  return length(max(d,0.0));
}
/*
float union(float d1, float d2) {
  return min(d1, d2);
}

float sub(float d1, float d2) {
  return max(-1.0 * d1, d2);
}

float intersect(float d1, float d2) {
  return max(d1, d2);
}
*/

float rayMarch(vec3 eye, vec3 dir) {
  float t = 0.01;
  int max_steps = 50;
  for (int i = 0; i < max_steps; i++) {
    vec3 p = eye + t * dir;
    float dist_sph = sdf_sphere(p, 5.0);
    //float dist_box = sdf_box(p, vec3(3.0,3.0,3.0));
    //float dist = union(dist_sph, dist_box);
    if (dist_sph < 0.001) {
      //at surface
      return t;
    }
    //move along ray
    t += dist_sph;
    if (t >= 1000.0) {
      //end
      return 1000.0;
    }
  }
  return 1000.0;
}

void main() {

  // convert to NDC screen coors
  vec4 s = vec4(-1.0 * (((gl_FragCoord.x / u_Dimensions.x) * 2.0) - 1.0),
                -1.0 * (1.0 - ((gl_FragCoord.y / u_Dimensions.y) * 2.0)), 1.0, 1.0);
  vec3 dir = rayCast(s);

  if (rayMarch(u_Eye, dir) < 1000.0) {
    //hit object
    out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
  }


  //out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
  
  //out_Col = vec4(0.5 * (fs_Pos + vec2(1.0)), 0.5 * (sin(u_Time * 3.14159 * 0.01) + 1.0), 1.0);

}
