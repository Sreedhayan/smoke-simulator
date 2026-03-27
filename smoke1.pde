// Hyper Realistic Smoke Simulation
// Fluid-based smoke using semi-Lagrangian advection

int N = 250;         // grid resolution (higher = more detail but slower)
float dt = 0.1;
float diffusion = 0.00001;
float viscosity = 0.00001;
float fadeSpeed = 0.003;

Fluid fluid;

void setup() {
  size(900, 900, P2D);
  fluid = new Fluid(N, diffusion, viscosity, dt);
  background(0);
}

void draw() {

  background(0);

  int cx = int(map(mouseX, 0, width, 10, N-10));
  int cy = int(map(mouseY, 0, height, 10, N-10));

  if (mousePressed) {
    for (int i = -3; i <= 3; i++) {
      for (int j = -3; j <= 3; j++) {

        fluid.addDensity(cx+i, cy+j, random(50,150));

        float angle = noise(frameCount*0.01, i, j) * TWO_PI;
        fluid.addVelocity(cx+i, cy+j,
                          cos(angle)*0.4,
                          sin(angle)*0.4);
      }
    }
  }

  fluid.step();
  fluid.render();
}

////////////////////////////////////////////////////////////

class Fluid {

  int size;
  float dt;
  float diff;
  float visc;

  float[] s;
  float[] density;

  float[] Vx;
  float[] Vy;

  float[] Vx0;
  float[] Vy0;

  Fluid(int size, float diffusion, float viscosity, float dt) {

    this.size = size;
    this.diff = diffusion;
    this.visc = viscosity;
    this.dt = dt;

    int total = (size+2)*(size+2);

    s = new float[total];
    density = new float[total];

    Vx = new float[total];
    Vy = new float[total];

    Vx0 = new float[total];
    Vy0 = new float[total];
  }

  int IX(int x, int y) {
    return x + (size+2) * y;
  }

  void addDensity(int x, int y, float amount) {
    density[IX(x,y)] += amount;
  }

  void addVelocity(int x, int y, float amountX, float amountY) {

    int index = IX(x,y);

    Vx[index] += amountX;
    Vy[index] += amountY;
  }

  void step() {

    diffuse(1, Vx0, Vx, visc);
    diffuse(2, Vy0, Vy, visc);

    project(Vx0, Vy0, Vx, Vy);

    advect(1, Vx, Vx0, Vx0, Vy0);
    advect(2, Vy, Vy0, Vx0, Vy0);

    project(Vx, Vy, Vx0, Vy0);

    diffuse(0, s, density, diff);
    advect(0, density, s, Vx, Vy);

    fadeDensity();
  }

  void fadeDensity() {
    for (int i=0;i<density.length;i++) {
      density[i] *= (1.0 - fadeSpeed);
    }
  }

  void render() {

    float cellW = width/(float)size;
    float cellH = height/(float)size;

    noStroke();

    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {

        float d = density[IX(i,j)];

        d = constrain(d,0,255);

        fill(d, d*0.95, d*0.9);

        rect(i*cellW, j*cellH, cellW+1, cellH+1);
      }
    }
  }

  ////////////////////////////////////////////////////////////

  void diffuse(int b, float[] x, float[] x0, float diff) {

    float a = dt * diff * size * size;

    for (int k=0;k<8;k++) {

      for (int i=1;i<=size;i++) {
        for (int j=1;j<=size;j++) {

          x[IX(i,j)] =
            (x0[IX(i,j)] +
             a * (
             x[IX(i-1,j)] +
             x[IX(i+1,j)] +
             x[IX(i,j-1)] +
             x[IX(i,j+1)]
             )) / (1 + 4*a);
        }
      }

      set_bnd(b, x);
    }
  }

  void advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY) {

    float dt0 = dt * size;

    for (int i=1;i<=size;i++) {
      for (int j=1;j<=size;j++) {

        float x = i - dt0 * velocX[IX(i,j)];
        float y = j - dt0 * velocY[IX(i,j)];

        if (x < 0.5) x = 0.5;
        if (x > size + 0.5) x = size + 0.5;
        int i0 = int(x);
        int i1 = i0 + 1;

        if (y < 0.5) y = 0.5;
        if (y > size + 0.5) y = size + 0.5;
        int j0 = int(y);
        int j1 = j0 + 1;

        float s1 = x - i0;
        float s0 = 1 - s1;

        float t1 = y - j0;
        float t0 = 1 - t1;

        d[IX(i,j)] =
          s0*(t0*d0[IX(i0,j0)] +
          t1*d0[IX(i0,j1)]) +
          s1*(t0*d0[IX(i1,j0)] +
          t1*d0[IX(i1,j1)]);
      }
    }

    set_bnd(b, d);
  }

  void project(float[] velocX, float[] velocY,
               float[] p, float[] div) {

    for (int i=1;i<=size;i++) {
      for (int j=1;j<=size;j++) {

        div[IX(i,j)] =
          -0.5*(
          velocX[IX(i+1,j)]
          -velocX[IX(i-1,j)]
          +velocY[IX(i,j+1)]
          -velocY[IX(i,j-1)]
          )/size;

        p[IX(i,j)] = 0;
      }
    }

    set_bnd(0, div);
    set_bnd(0, p);

    for (int k=0;k<20;k++) {

      for (int i=1;i<=size;i++) {
        for (int j=1;j<=size;j++) {

          p[IX(i,j)] =
            (div[IX(i,j)]
            + p[IX(i-1,j)]
            + p[IX(i+1,j)]
            + p[IX(i,j-1)]
            + p[IX(i,j+1)]) / 4;
        }
      }

      set_bnd(0,p);
    }

    for (int i=1;i<=size;i++) {
      for (int j=1;j<=size;j++) {

        velocX[IX(i,j)] -=
          0.5*(p[IX(i+1,j)]
          -p[IX(i-1,j)])*size;

        velocY[IX(i,j)] -=
          0.5*(p[IX(i,j+1)]
          -p[IX(i,j-1)])*size;
      }
    }

    set_bnd(1, velocX);
    set_bnd(2, velocY);
  }

  void set_bnd(int b, float[] x) {

    for (int i=1;i<=size;i++) {

      x[IX(0,i)] =
        b==1 ? -x[IX(1,i)] : x[IX(1,i)];

      x[IX(size+1,i)] =
        b==1 ? -x[IX(size,i)] : x[IX(size,i)];

      x[IX(i,0)] =
        b==2 ? -x[IX(i,1)] : x[IX(i,1)];

      x[IX(i,size+1)] =
        b==2 ? -x[IX(i,size)] : x[IX(i,size)];
    }

    x[IX(0,0)] =
      0.5*(x[IX(1,0)] + x[IX(0,1)]);

    x[IX(0,size+1)] =
      0.5*(x[IX(1,size+1)] + x[IX(0,size)]);

    x[IX(size+1,0)] =
      0.5*(x[IX(size,0)] + x[IX(size+1,1)]);

    x[IX(size+1,size+1)] =
      0.5*(x[IX(size,size+1)] + x[IX(size+1,size)]);
  }
}
