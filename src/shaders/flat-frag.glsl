#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

vec3 light;

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

    //set light vector for lambert shading
    light = vec3((inverse(view) * vec4(0.0, 0.0, 0.0, 1.0)) - vec4(fs_Pos, 1.0, 1.0));

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

float sdf_torus(vec3 p, vec2 t) {
  vec2 q = vec2(length(p.xz) - t.x, p.y);
  return length(q) - t.y;
}

float union_op(float d1, float d2) {
  return min(d1, d2);
}

float sub_op(float d1, float d2) {
  return max(-1.0 * d1, d2);
}

float sect_op(float d1, float d2) {
  return max(d1, d2);
}

mat4 rot_matrix(vec3 r) {
  //convert angle to radians
  r = (r * 3.14159) / 180.0;

  mat4 m;
  m[0][0] = 1.0;
  m[1][1] = 1.0;
  m[2][2] = 1.0;
  m[3][3] = 1.0;
  mat4 rot_x = m;
  rot_x[1][1] = cos(r.x);
  rot_x[1][2] = sin(r.x);
  rot_x[2][1] = -1.0 * sin(r.x);
  rot_x[2][2] = cos(r.x);
  mat4 rot_y = m;
  rot_y[0][0] = cos(r.y);
  rot_y[2][2] = cos(r.y);
  rot_y[2][0] = sin(r.y);
  rot_y[0][2] = -1.0 * sin(r.y);
  mat4 rot_z = m;
  rot_z[0][0] = cos(r.z);
  rot_z[1][1] = cos(r.z);
  rot_z[1][0] = -1.0 * sin(r.z);
  rot_z[0][1] = sin(r.z);

  return rot_x * rot_y * rot_z;
  //return rot_x;
  //return m;
}

vec3 trans_pt(vec3 p, vec3 t) {
  return p + t;
}

vec3 rot_op(vec3 r, vec3 p) {
  mat4 tmat = rot_matrix(-1.0 * r);
  vec3 p_rot = vec3(tmat * vec4(p, 1.0));
  return p_rot;
}

float sceneSDF(vec3 p) {
  //apply rotation + translation to box
  vec3 rot_box = vec3(0.0, 0.0, 60.0);
  vec3 p_box = rot_op(rot_box, p);
  vec3 t_box = vec3(2.3, 0.0, 0.0);

  //apply scale + translation to sphere
  vec3 t_sph = vec3(0.0, 0.5, 0.0);
  float dist_sph = sdf_sphere(trans_pt(p, t_sph) * 1.0, 3.0) / 1.0;

  float dist_box = sdf_box(trans_pt(p_box, t_box) * 0.9, vec3(3.0, 3.0, 3.0)) / 0.9;
  float dist_cup = sect_op(dist_box, dist_sph);

  float dist_sph2 = sdf_sphere(p * 1.3, 5.0) / 1.3;
  vec3 t_box2 = vec3(0.0, 4.0, 0.0);
  float dist_box2 = sdf_box(trans_pt(p, t_box2) * 0.6, vec3(3.0, 3.0, 3.0)) / 0.6;
  float dist_cup2 = sect_op(dist_box2, dist_sph2);

  //cup handle
  vec3 rot_torus = vec3(90.0, 0.0, 0.0);
  vec3 p_torus = rot_op(rot_torus, p);
  vec3 t_tor = vec3(2.5, 0.0, 0.8);
  float dist_torus = sdf_torus(trans_pt(p_torus, t_tor), vec2(1.2, 0.2));

  dist_cup = sect_op(dist_box, dist_sph);
  float dist = union_op(dist_cup, dist_torus);
  return dist;
}

vec3 estNormal(vec3 p) {
  float eps = 0.001;
  vec3 nor = vec3(sceneSDF(vec3(p.x + eps, p.y, p.z)) - sceneSDF(vec3(p.x - eps, p.y, p.z)),
                  sceneSDF(vec3(p.x, p.y + eps, p.z)) - sceneSDF(vec3(p.x, p.y - eps, p.z)),
                  sceneSDF(vec3(p.x, p.y, p.z + eps)) - sceneSDF(vec3(p.x, p.y, p.z - eps)));
  return normalize(nor);
}

vec4 rayMarch(vec3 eye, vec3 dir) { 
  //rayMarch returns (t,vec3 nor)
  float t = 0.01;
  int max_steps = 50;
  for (int i = 0; i < max_steps; i++) {
    vec3 p = eye + t * dir;

    float dist = sceneSDF(p);
    vec3 nor = estNormal(p);

    if (dist < 0.001) {
      //at surface
      return vec4(t, nor);
    }
    //move along ray
    t += dist;
    if (t >= 1000.0) {
      //end
      return vec4(t, 0.0, 0.0, 0.0);
    }
  }
  t = 1000.0;
  return vec4(t, 0.0, 0.0, 0.0);
}

void main() {

  // convert to NDC screen coors
  vec4 s = vec4(-1.0 * (((gl_FragCoord.x / u_Dimensions.x) * 2.0) - 1.0),
                -1.0 * (1.0 - ((gl_FragCoord.y / u_Dimensions.y) * 2.0)), 1.0, 1.0);
  vec3 dir = rayCast(s);

  vec4 march = rayMarch(u_Eye, dir);
  if (march[0] < 1000.0) {
    //hit object
    out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
    
    vec4 diffuseColor = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
    
    // Calculate the diffuse term for Lambert shading
    vec3 nor = vec3(march[1],march[2],march[3]);
    
    float diffuseTerm = dot(normalize(nor), normalize(light));
    // Avoid negative lighting values
    diffuseTerm = clamp(diffuseTerm, 0.0, 1.0);
    float ambientTerm = 0.2;
    float lightIntensity = diffuseTerm + ambientTerm;
    // Compute final shaded color
    out_Col = vec4(diffuseColor.rgb * lightIntensity, 1.0);
  }


  //out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
  
  //out_Col = vec4(0.5 * (fs_Pos + vec2(1.0)), 0.5 * (sin(u_Time * 3.14159 * 0.01) + 1.0), 1.0);

}
