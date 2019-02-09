#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;
uniform vec4 u_Color;
uniform float u_Anim;

in vec2 fs_Pos;
out vec4 out_Col;

vec3 light;
vec3 nor;
float worl;

vec3 anim_trans;
vec3 anim_angle;

vec3 rayCast(vec4 s) {
    // multiply by far clip plane
    float far_clip = 1000.0;
    float near_clip = 0.1;
    s *= far_clip;

    // multiply by inverse projection matrix
    mat4 proj;
    float fov = 45.0;
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

    // set up ray
    vec3 origin = u_Eye;
    vec3 dir = normalize(vec3(s) - u_Eye);

    // set light vector for shading
    light = vec3((inverse(view) * vec4(0.0, 0.0, 0.0, 1.0)) - vec4(fs_Pos, 1.0, 1.0));

    return dir;
}

vec2 random2(vec2 p, vec2 seed) {
  return fract(sin(vec2(dot(p + seed, vec2(311.7, 127.1)), dot(p + seed, vec2(269.5, 183.3)))) * 85734.3545);
}

float WorleyNoise(vec2 uv) {
    // Tile the space
    vec2 uvInt = floor(uv);
    vec2 uvFract = fract(uv);

    float minDist = 1.0; // Minimum distance initialized to max.

    // Search all neighboring cells and this cell for their point
    for(int y = -1; y <= 1; y++) {
        for(int x = -1; x <= 1; x++) {
            vec2 neighbor = vec2(float(x), float(y));

            // Random point inside current neighboring cell
            vec2 point = random2(uvInt + neighbor, vec2(10.0));

            // Compute the distance b/t the point and the fragment
            // Store the min dist thus far
            vec2 diff = neighbor + point - uvFract;
            float dist = length(diff);
            minDist = min(minDist, dist);
        }
    }
    return minDist;
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

float sdf_cylin( vec3 p, vec2 h) {
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
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
  // convert angle to radians
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
}

vec3 trans_pt(vec3 p, vec3 t) {
  return p + t;
}

vec3 rot_op(vec3 r, vec3 p) {
  mat4 tmat = rot_matrix(-1.0 * r);
  vec3 p_rot = vec3(tmat * vec4(p, 1.0));
  return p_rot;
}

float cupSDF(vec3 p) {
  // apply rotation + translation to box point
  vec3 rot_box = vec3(0.0, 55.0, 60.0);
  vec3 t_box = vec3(0.5, 0.0, 0.0);
  t_box += anim_trans;
  vec3 p_box = trans_pt(p, t_box);

  // apply translation to sphere
  vec3 t_sph = vec3(0.0, 0.5, 0.0);
  t_sph += anim_trans;
  float dist_sph = sdf_sphere(trans_pt(p, t_sph), 3.0);

  // create basic cup shape with sphere + box
  float dist_box = sdf_box(rot_op(rot_box, p_box) * 0.8, vec3(1.2, 3.0, 3.0)) / 0.8;
  float dist_cup = sect_op(dist_box, dist_sph);

  // create cup handle with torus
  vec3 rot_torus = vec3(65.0, 0.0, 0.0);
  vec3 t_tor = vec3(2.5, 0.0, -0.25);
  t_tor += anim_trans;
  vec3 p_torus = trans_pt(p, t_tor);
  float dist_torus = sdf_torus(rot_op(rot_torus, p_torus) * 0.75, vec2(1.2, 0.2)) / 0.75;

  // carve out inside of cup
  vec3 t_sph2 = vec3(0.0, 0.5, 0.0);
  t_sph2 += anim_trans;
  float dist_sph2 = sdf_sphere(trans_pt(p, t_sph2) * 1.0, 2.8) / 1.0;
  vec3 t_box2 = vec3(0.0, 0.0, 5.0);
  t_box2 += anim_trans;
  vec3 p_box2 = trans_pt(p, t_box2);
  float dist_box2 = sdf_box(rot_op(rot_box, p_box2) * 0.3, vec3(1.2, 3.0, 3.0)) / 0.3;
  vec3 t_box3 = vec3(0.0, 1.7, -0.8);
  t_box3 += anim_trans;
  vec3 p_box3 = trans_pt(p, t_box3);
  float dist_box3 = sdf_box(rot_op(rot_box, p_box3) * 0.7, vec3(0.2, 2.5, 2.5)) / 0.7;

  dist_cup = sect_op(dist_box, dist_sph);
  // connect cup to handle
  float dist = union_op(dist_cup, dist_torus);
  // carve out inside of cup
  float dist_inside = sect_op(dist_sph2, dist_box2);
  dist = sub_op(dist_sph2, dist);
  // add in bottom of cup
  dist_inside = sect_op(dist_sph2, dist_box3);
  dist = union_op(dist_inside, dist);

  return dist;
}

float coffeeSDF(vec3 p) {
  // SDF for the liquid inside cup
  // sphere
  vec3 t_sph2 = vec3(0.0, 0.5, 0.0);
  t_sph2 += anim_trans;
  float dist_sph2 = sdf_sphere(trans_pt(p, t_sph2) * 1.0, 2.8) / 1.0;
  // box
  vec3 t_box3 = vec3(0.0, 4.8, 0.0);
  t_box3 += anim_trans;
  float dist_box3 = sdf_box(trans_pt(p, t_box3) * 1.0, vec3(3.0, 3.0, 3.0)) / 1.0;
  float dist = sect_op(dist_sph2, dist_box3);

  // apply worley noise for water effect
  worl = WorleyNoise(2.0 * p.xz);

  return dist;
}

float charSDF(vec3 p) {
  // char head
  vec3 t_sph2 = vec3(-0.5, 1.8, 0.0);
  t_sph2 += anim_trans;
  vec3 p_sph2 = trans_pt(p, t_sph2);
  float dist_head = sdf_sphere(p_sph2 * 1.0, 0.4) / 1.0;

  // char arms
  vec3 t_arm1 = vec3(0.1, 2.0, 0.0);
  vec3 r_arm1 = vec3(0.0, 0.0, 35.0);
  t_arm1 += anim_trans;
  vec3 p_arm1 = trans_pt(p, t_arm1);
  r_arm1 += anim_angle;
  float dist_arm1 = sdf_cylin(rot_op(r_arm1, p_arm1), vec2(0.08, 0.7));
  vec3 t_arm2 = vec3(-1.1, 2.0, 0.0);
  vec3 r_arm2 = vec3(0.0, 0.0, -35.0);
  t_arm2 += anim_trans;
  vec3 p_arm2 = trans_pt(p, t_arm2);
  r_arm2 += anim_angle;
  float dist_arm2 = sdf_cylin(rot_op(r_arm2, p_arm2), vec2(0.08, 0.7));

  // connect head and arms into same SDF
  float dist = union_op(dist_head, dist_arm1);
  dist = union_op(dist, dist_arm2);

  return dist;
}

float charEyeSDF(vec3 p) { 
  // char eyes 
  vec3 t_sph1 = vec3(-0.3, 1.5, 0.3);
  t_sph1 += anim_trans;
  vec3 p_sph1 = trans_pt(p, t_sph1);
  float dist_eye1 = sdf_sphere(p_sph1, 0.05);
  vec3 t_sph2 = vec3(-0.7, 1.5, 0.3);
  t_sph2 += anim_trans;
  float dist_eye2 = sdf_sphere(trans_pt(p, t_sph2), 0.05);

  //char mouth
  vec3 t_tor = vec3(-0.5, 1.5, 0.5);
  t_tor += anim_trans;
  vec3 p_tor = trans_pt(p, t_tor);
  vec3 r_tor = vec3(65.0, 0.0, 0.0);
  float dist_tor = sdf_torus(rot_op(r_tor, p_tor), vec2(0.08, 0.02));
  vec3 t_sph3 = vec3(-0.5, 1.65, 0.3);
  t_sph3 += anim_trans;
  float dist_sph3 = sdf_sphere(trans_pt(p, t_sph3), 0.25);

  // connect eyes to same SDF
  float dist = union_op(dist_eye1, dist_eye2);
  // create mouth by taking part of torus
  float dist_mouth = sect_op(dist_sph3, dist_tor);
  // combine eyes and mouth to same SDF
  dist = union_op(dist, dist_mouth);
  return dist;
}

vec3 estNormalCup(vec3 p) {
  // find normal of cup points
  float eps = 0.001;
  vec3 nor_c = vec3(cupSDF(vec3(p.x + eps, p.y, p.z)) - cupSDF(vec3(p.x - eps, p.y, p.z)),
                  cupSDF(vec3(p.x, p.y + eps, p.z)) - cupSDF(vec3(p.x, p.y - eps, p.z)),
                  cupSDF(vec3(p.x, p.y, p.z + eps)) - cupSDF(vec3(p.x, p.y, p.z - eps)));
  return normalize(nor_c);
}

vec3 estNormalCoffee(vec3 p) {
  // find normal of liquid points
  float eps = 0.001;
  vec3 nor_co = vec3(coffeeSDF(vec3(p.x + eps, p.y, p.z)) - coffeeSDF(vec3(p.x - eps, p.y, p.z)),
                  coffeeSDF(vec3(p.x, p.y + eps, p.z)) - coffeeSDF(vec3(p.x, p.y - eps, p.z)),
                  coffeeSDF(vec3(p.x, p.y, p.z + eps)) - coffeeSDF(vec3(p.x, p.y, p.z - eps)));
  return normalize(nor_co);
}

vec2 rayMarch(vec3 eye, vec3 dir) { 
  // rayMarch returns (t, object id)
  float t = 0.01;
  int max_steps = 1000;
  vec3 p = eye + t * dir;
  for (int i = 0; i < max_steps; i++) {
    p = eye + t * dir;

    float dist = cupSDF(p);
    float dist2 = coffeeSDF(p);
    float dist3 = charSDF(p);
    float dist4 = charEyeSDF(p);

    nor = estNormalCup(p);

    if (dist < 0.00001) {
      // at cup surface
      nor = estNormalCup(p);
      return vec2(t, 1.0);  
      // move along ray
      t += dist;
    }
    else if (dist2 < 0.00001) {
      // at coffee surface
      nor = estNormalCoffee(p);
      return vec2(t, 2.0);
      // move along ray
      t += dist2;
    }
    else if (dist3 < 0.00001) {
      // at character surface
      // flat shading; no nor set
      return vec2(t, 3.0);
      // move along ray
      t += dist3;
    }
    else if (dist4 < 0.00001) {
      // at character eye surface
      // flat shading; no nor set
      return vec2(t, 4.0);
      // move alone ray
      t += dist4;
    }
    else {
      // increment by smallest distance
      float dist_min = min(dist, dist2);
      dist_min = min(dist_min, dist3);
      dist_min = min(dist_min, dist4);
      t += dist_min;
    }
  
    if (t >= 1000.0) {
      // end
      return vec2(t, 0.0);
    }
  }
  t = 1000.0;
  return vec2(t, 0.0);
}

bool rayBoxIntersection(vec3 origin, vec3 dir, vec3 min, vec3 max) {
  // check intersection of ray with cube for bounding box purposes
  float near = -1.0 * (1.0 / 0.0);
  float far = (1.0 / 0.0);
  float t0;
  float t1;
  for (int i = 0; i < 3; i++) {
    if (dir[i] == 0.0) {
      if (origin[i] < min[i] || origin[i] > max[i]) {
        return false;
      }
    }
    t0 = (min[i] - origin[i]) / dir[i];
    t1 = (max[i] - origin[i]) / dir[i];
    if (t0 > t1) {
      float temp = t0;
      t0 = t1;
      t1 = temp;
    }
    if (t0 > near) {
      near = t0;
    }
    if (t1 < far) {
      far = t1;
    }
  }
  if (near > far) {
    return false;
  }
  else {
    return true;
  }
}

float ease_in_out_quadratic(float t) {
  if (t < 0.5) {
    return (t * t * 4.0) / 2.0;
  }
  else {
    return 1.0 - (4.0 - (8.0 * t) + (t * t)) / 2.0;
  }
}

float ease_linear(float t, float b, float c, float d) {
  return c * (t / d) + b;
}

vec3 colorConvert(vec3 rgb) {
  vec3 new_color;
  new_color[0] = (rgb[0] / 255.0);
  new_color[1] = (rgb[1] / 255.0);
  new_color[2] = (rgb[2] / 255.0);
  return new_color;
}

void main() {
  // convert to NDC screen coors
  vec4 s = vec4(-1.0 * (((gl_FragCoord.x / u_Dimensions.x) * 2.0) - 1.0),
                -1.0 * (1.0 - ((gl_FragCoord.y / u_Dimensions.y) * 2.0)), 1.0, 1.0);
  vec3 dir = rayCast(s);

  // set up global animation values
  float time = float(u_Time);
  anim_trans = vec3(0.0, 0.5 * ease_in_out_quadratic(sin(time * 3.14159 * 0.01)), 0.0);
  // character animation dependent on passed in value
  if (u_Anim == 0.0) {
    // no animation
    anim_angle = vec3(0.0);
  }
  else {
    // animate character arms
    anim_angle = vec3(0.0, 0.0, (ease_linear(time, 0.0, 5.0, 3.0)));
    anim_angle[2] = mod(anim_angle[2], 120.0);
    if (anim_angle[2] > 60.0) {
    anim_angle[2] = 120.0 - anim_angle[2];
    }
  }
  
  // set up bounding box hierarchy
  // bounds for main sphere for cup
  vec3 center = vec3(0.0,0.0,0.0);
  vec3 box_cup_sph_min = center - vec3(0.0, 0.5, 0.0) - vec3(3.0);
  vec3 box_cup_sph_max = center - vec3(0.0, 0.5, 0.0) + vec3(3.0);
  box_cup_sph_min -= anim_trans;
  box_cup_sph_max -= anim_trans;

  // bounds for torus handle
  vec3 box_cup_tor_min = center - vec3(2.5, 0.0, -0.25) - vec3(1.2);
  box_cup_tor_min *= (1.0 / 0.75);
  vec3 box_cup_tor_max = center - vec3(2.5, 0.0, -0.25) + vec3(1.2);
  box_cup_tor_max *= (1.0 / 0.75);
  box_cup_tor_min -= anim_trans;
  box_cup_tor_max -= anim_trans;

  // root contains both cup box and handle box
  vec3 box_min = min(box_cup_sph_min, box_cup_tor_min);
  vec3 box_max = max(box_cup_sph_max, box_cup_tor_max);

  bool bound_test = rayBoxIntersection(u_Eye, dir, box_min, box_max);
  if (bound_test) {
    // in root bounding box
    bool bound_test_cup = rayBoxIntersection(u_Eye, dir, box_cup_sph_min, box_cup_sph_max);
    bool bound_test_tor = rayBoxIntersection(u_Eye, dir, box_cup_tor_min, box_cup_tor_max);
    if (bound_test_cup || bound_test_tor) {
      // in cup or handle bounding boxes
      // rayMarch test
      vec2 march = rayMarch(u_Eye, dir);

      if (march[0] < 1000.0) {
        // Hit object
        vec4 diffuseColor;
        if (march[1] == 1.0) {
          // Hit cup
          diffuseColor = vec4(colorConvert(vec3(194.0, 187.0, 240.0)), 1.0);
        }
        else if (march[1] == 2.0) {
          // Hit coffee
          diffuseColor = vec4(worl);
          diffuseColor += vec4(colorConvert(vec3(0.0, 100.0, 200.0)), 1.0);
          while (diffuseColor[2] < 235.0) {
            diffuseColor += vec4(colorConvert(vec3(0.0, 100.0, 200.0)), 1.0);
          }
        }
        else if (march[1] == 3.0) {
          // Hit char
          diffuseColor = vec4(u_Color / 255.0);
        }
        else if (march[1] == 4.0) {
          // Hit char features
          diffuseColor = vec4(colorConvert(vec3(255.0, 223.0, 109.0)), 1.0);
        }

        // Calculate diffuse term for shading
        float diffuseTerm = dot(normalize(nor), normalize(light));
        // Avoid negative lighting values
        diffuseTerm = clamp(diffuseTerm, 0.0, 1.0);
        float ambientTerm = 0.2;
        float lightIntensity = diffuseTerm + ambientTerm;

        // Implement specular light
        vec4 H;
        for (int i = 0; i < 4; i++) {
            H[i] = (light[i] + u_Eye[i]) / 2.0;
        }
        float specularIntensity = max(pow(dot(normalize(H), normalize(vec4(nor,1.0))), 1.5), 0.0);

        // Compute final shaded color
        out_Col = vec4(diffuseColor.rgb * (lightIntensity + specularIntensity), 1.0);
      }
      else {
        // bg color
        out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
      }
    }
    else {
      // bg color
      out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
    }
  }
  else {
    // bg color
    out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
  }
}
