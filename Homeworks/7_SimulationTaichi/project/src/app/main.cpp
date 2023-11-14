//  88-Line 2D Moving Least Squares Material Point Method (MLS-MPM)
// [Explained Version by David Medina]

// Uncomment this line for image exporting functionality
#define TC_IMAGE_IO

// Note: You DO NOT have to install taichi or taichi_mpm.
// You only need [taichi.h] - see below for instructions.
#include "taichi.h"

// ---------- Add random ----------
#include <ctime>
#include <cstdlib>
// ---------- End ----------

using namespace taichi;
using Vec = Vector2;
using Mat = Matrix2;

// Window
const int window_size = 800;

// Grid resolution (cells)
const int n = 80;

const real dt = 1e-4_f;
const real frame_dt = 1e-3_f;
const real dx = 1.0_f / n;
const real inv_dx = 1.0_f / dx;

// Snow material properties
const auto particle_mass = 1.0_f;
const auto vol = 1.0_f;        // Particle Volume
const auto hardening = 10.0_f; // Snow hardening factor
const auto E = 1e4_f;          // Young's Modulus
const auto nu = 0.2_f;         // Poisson ratio
const bool plastic = true;

// Initial Lamé parameters
const real mu_0 = E / (2 * (1 + nu));
const real lambda_0 = E * nu / ((1 + nu) * (1 - 2 * nu));

struct Particle {
	// Position and velocity
	Vec x, v;
	// Deformation gradient
	Mat F;
	// Affine momentum from APIC
	Mat C;
	// Determinant of the deformation gradient (i.e. volume)
	real Jp;
	// Color
	int c;
	//---------- Add Particle Type ----------
	// Type(0:Snow,1:Liquid,2:Jelly)
	int type_p;

	Particle(Vec x, int c, int type_p = 0, Vec v = Vec(0)) :
		x(x),
		v(v),
		F(1),
		C(0),
		Jp(1),
		c(c),
		type_p(type_p) {}
	//---------- End ----------
};

std::vector<Particle> particles;

// Vector3: [velocity_x, velocity_y, mass]
Vector3 grid[n + 1][n + 1];

void advance(real dt) {
	// Reset grid
	std::memset(grid, 0, sizeof(grid));

	// P2G
	for (auto& p : particles) {
		// element-wise floor
		Vector2i base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();

		Vec fx = p.x * inv_dx - base_coord.cast<real>();

		// Quadratic kernels [http://mpm.graphics Eqn. 123, with x=fx, fx-1,fx-2]
		Vec w[3] = {
		  Vec(0.5) * sqr(Vec(1.5) - fx),
		  Vec(0.75) - sqr(fx - Vec(1.0)),
		  Vec(0.5) * sqr(fx - Vec(0.5))
		};

		// ---------- Add new particle types ----------
		// Compute current Lamé parameters [http://mpm.graphics Eqn. 86]
		auto e = std::exp(hardening * (1.0f - p.Jp));
		if (p.type_p == 2) e = 1;     // Jelly
		auto mu = mu_0 * e;
		auto lambda = lambda_0 * e;
		if (p.type_p == 1) mu = 0;      // Liquid
		// ---------- End ----------

		// Current volume
		real J = determinant(p.F);

		// Polar decomposition for fixed corotated model
		Mat r, s;
		polar_decomp(p.F, r, s);

		// [http://mpm.graphics Paragraph after Eqn. 176]
		real Dinv = 4 * inv_dx * inv_dx;
		// [http://mpm.graphics Eqn. 52]
		auto PF = (2 * mu * (p.F - r) * transposed(p.F) + lambda * (J - 1) * J);

		// Cauchy stress times dt and inv_dx
		auto stress = -(dt * vol) * (Dinv * PF);

		// Fused APIC momentum + MLS-MPM stress contribution
		// See http://taichi.graphics/wp-content/uploads/2019/03/mls-mpm-cpic.pdf
		// Eqn 29
		auto affine = stress + particle_mass * p.C;

		// P2G
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				auto dpos = (Vec(i, j) - fx) * dx;
				// Translational momentum
				Vector3 mass_x_velocity(p.v * particle_mass, particle_mass);
				grid[base_coord.x + i][base_coord.y + j] += (
					w[i].x * w[j].y * (mass_x_velocity + Vector3(affine * dpos, 0))
					);
			}
		}
	}

	// For all grid nodes
	real mu_v = 0.5;
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			auto& g = grid[i][j];
			// No need for epsilon here
			if (g[2] > 0) {
				// Normalize by mass
				g /= g[2];
				// Gravity
				// ---------- For Tetris ----------
				g += dt * Vector3(0, -200, 0);
				if (g[1] > 0) g[1] *= exp(-1 / mu_v * (j + 1));
				// ---------- End ----------

				// ---------- For Planet ----------
				//auto r = Vec(n / 2.0_f - i, n / 2.0_f - j);
				//auto f = real(5000 / pow(r.length() > 0.1 ? r.length() : 0.1, 3)) * r;
				//g += dt * Vector3(f.x, f.y, 0);
				// ---------- End ----------

				// boundary thickness
				real boundary = 0.05;
				// Node coordinates
				real x = (real)i / n;
				real y = real(j) / n;

				// Sticky boundary
				if (x < boundary || x > 1 - boundary || y > 1 - boundary) {
					g = Vector3(0);
				}
				// Separate boundary
				if (y < boundary) {
					g[1] = std::max(0.0f, g[1]);
					// ---------- For Tetris ----------
					g[0] *= mu_v;
					// ---------- End ----------
				}
			}
		}
	}

	// G2P
	for (auto& p : particles) {
		// element-wise floor
		Vector2i base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();
		Vec fx = p.x * inv_dx - base_coord.cast<real>();
		Vec w[3] = {
					Vec(0.5) * sqr(Vec(1.5) - fx),
					Vec(0.75) - sqr(fx - Vec(1.0)),
					Vec(0.5) * sqr(fx - Vec(0.5))
		};

		p.C = Mat(0);
		p.v = Vec(0);

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				auto dpos = (Vec(i, j) - fx);
				auto grid_v = Vec(grid[base_coord.x + i][base_coord.y + j]);
				auto weight = w[i].x * w[j].y;
				// Velocity
				p.v += weight * grid_v;
				// APIC C
				p.C += 4 * inv_dx * Mat::outer_product(weight * grid_v, dpos);
			}
		}

		// Advection
		p.x += dt * p.v;

		// MLS-MPM F-update
		auto F = (Mat(1) + dt * p.C) * p.F;

		// ---------- Add new particle types ----------
		if (p.type_p == 0)
		{
			Mat svd_u, sig, svd_v;
			svd(F, svd_u, sig, svd_v);

			// Snow Plasticity
			for (int i = 0; i < 2 * int(plastic); i++) {
				sig[i][i] = clamp(sig[i][i], 1.0f - 2.5e-2f, 1.0f + 7.5e-3f);
			}

			real oldJ = determinant(F);
			F = svd_u * sig * transposed(svd_v);

			real Jp_new = clamp(p.Jp * oldJ / determinant(F), 0.6f, 20.0f);

			p.Jp = Jp_new;
			p.F = F;
		}
		else if (p.type_p == 1)
		{
			// Liquid
			p.F = Mat(1) * sqrt(determinant(F));
		}
		else if (p.type_p == 2)
		{
			// Jelly
			p.F = F;
		}
		// ---------- End ----------
	}
}

// Seed particles with position and color
void add_object(Vec center, int c, int num, int type_p = 0) {
	// Randomly sample 100 particles in the square
	for (int i = 0; i < num; i++) {
		particles.push_back(Particle((Vec::rand() * 2.0f - Vec(1)) * 0.03f + center, c, type_p));
	}
}

//---------- Add Tetris ----------
void add_object_tetris(Vec center, int c, int num, int type_tetris = 0, int type_p = 0)
{
	add_object(center, c, num, type_p);
	add_object(center + Vec(0, -0.06), c, num, type_p);

	switch (type_tetris)
	{
	case 0:
		add_object(center + Vec(0.06, -0.06), c, num, type_p);
		add_object(center + Vec(-0.06, -0.06), c, num, type_p);
		break;
	case 1:
		add_object(center + Vec(0.06, -0.06), c, num, type_p);
		add_object(center + Vec(-0.06, 0), c, num, type_p);
		break;
	case 2:
		add_object(center + Vec(0.06, 0), c, num, type_p);
		add_object(center + Vec(-0.06, -0.06), c, num, type_p);
		break;
	case 3:
		add_object(center + Vec(0, -0.12), c, num, type_p);
		add_object(center + Vec(0, -0.18), c, num, type_p);
		break;
	case 4:
		add_object(center + Vec(0.06, -0.06), c, num, type_p);
		add_object(center + Vec(0.06, 0), c, num, type_p);
		break;
	case 5:
		add_object(center + Vec(0.06, -0.06), c, num, type_p);
		add_object(center + Vec(0.12, -0.06), c, num, type_p);
		break;
	case 6:
		add_object(center + Vec(-0.06, -0.06), c, num, type_p);
		add_object(center + Vec(-0.12, -0.06), c, num, type_p);
		break;
	default:
		break;
	}
}

void add_tetris(int num, int type_p = -1, int type_tetris = -1, Vec center = Vec(-1))
{
	if (type_p == -1)
		type_p = (int)(std::rand() % 3);
	if (type_tetris == -1)
		type_tetris = (int)(std::rand() % 7);
	if (center.x == -1)
		center = Vec(std::rand() % 10 * 0.06 + 0.23, 0.90);
	switch (type_tetris)
	{
	case 0:
		add_object_tetris(center, 0x6C3365, num, type_tetris, type_p);
		break;
	case 1:
		add_object_tetris(center, 0xFF0000, num, type_tetris, type_p);
		break;
	case 2:
		add_object_tetris(center, 0x28FF28, num, type_tetris, type_p);
		break;
	case 3:
		add_object_tetris(center, 0x00FFFF, num, type_tetris, type_p);
		break;
	case 4:
		add_object_tetris(center, 0xFFFF37, num, type_tetris, type_p);
		break;
	case 5:
		add_object_tetris(center, 0x0000E3, num, type_tetris, type_p);
		break;
	case 6:
		add_object_tetris(center, 0xFFD306, num, type_tetris, type_p);
		break;
	default:
		break;
	}
}
//---------- End ----------

// ---------- Add Clear ----------
void tetris_clear(int num)
{
	int count = 0;
	for (auto& p : particles)
	{
		if (p.x.y < 0.11)
			count++;
	}
	if (count > 0.8 * num / (0.06 * 0.06) * 0.9 * 0.06)
	{
		for (auto iter = particles.begin(); iter < particles.end();)
		{
			if (iter->x.y < 0.11)
			{
				particles.erase(iter);
			}
			else
				iter++;
		}
	}
}
// ---------- End ----------

// ---------- Add Planet ----------
void add_planet(Vec center, int c, int type_p, real omega = 0.3, int num = 2000, real r = 0.1, Vec v = Vec(0))
{
	//for (int i = 0; i < 2000; i++) {
	//    auto pos = Vec::rand() * 2.0f - Vec(1);
	//    particles.push_back(Particle(pos * 0.1f + center, c, type_p, omega * Vec(pos.y, -pos.x)));
	//}
	int i = 0;
	while (i < num) {
		auto pos = (Vec::rand() * 2.0_f - Vec(1)) * r;
		if (pos.x * pos.x + pos.y * pos.y < r * r) {
			particles.push_back(Particle(pos + center, c, type_p, omega * Vec(pos.y, -pos.x) + v));
			i++;
		}
	}
}
// ---------- End ----------

void sample(int count)
{
	switch (count)
	{
	case 0:
		add_object_tetris(Vec(0.29, 0.80), 0xFF0000, 200, 1, 0);
		break;
	case 1:
		add_object_tetris(Vec(0.59, 0.80), 0xFF0000, 200, 1, 2);
		break;
	case 2:
		//add_object_tetris(Vec(0.69, 0.80), 0x00FFFF, 200, 2, 1);
		break;
	default:
		break;
	}
}
int main() {
	GUI gui("Real-time 2D MLS-MPM", window_size, window_size);
	auto& canvas = gui.get_canvas();

	srand((int)time(0));
	int num = 100;
	int frame = 0;
	int count = 0;

	// ---------- For Planet ----------
	//add_planet(Vec(0.50, 0.50), 0xFFFF37, 1, 60);
	//add_planet(Vec(0.75, 0.50), 0x00FFFF, 2, 0.3, 200, 0.03, Vec(0, -1.5));
	//add_planet(Vec(0.85, 0.50), 0xFF0000, 0, 0.2, 500, 0.03, Vec(0, -1.5));
	// ---------- End ----------

	// Main Loop
	for (int step = 0;; step++) {
		// Advance simulation
		advance(dt);

		// Visualize frame
		if (step % int(frame_dt / dt) == 0) {
			// Clear background
			canvas.clear(0x112F41);
			// Box
			canvas.rect(Vec(0.04), Vec(0.96)).radius(2).color(0x4FB99F).close();
			// Particles
			for (auto p : particles) {
				canvas.circle(p.x).radius(2).color(p.c);
			}
			// Update image
			gui.update();

			// ---------- For Tetris ----------
			if (step % int(150 * frame_dt / dt) == 0 && count < 10) {
				tetris_clear(num);
				//sample(count);
				add_tetris(num);
				//add_object_tetris(Vec(0.50, 0.85), 0xFFFF37, num, 3, 0);
				count++;
			}
			// ---------- End ----------
			// Write to disk (optional)
			//canvas.img.write_as_image(fmt::format("tmp/tetris2/{:05d}.png", frame++));
		}
	}
}

/* -----------------------------------------------------------------------------
** Reference: A Moving Least Squares Material Point Method with Displacement
			  Discontinuity and Two-Way Rigid Body Coupling (SIGGRAPH 2018)

  By Yuanming Hu (who also wrote this 88-line version), Yu Fang, Ziheng Ge,
		   Ziyin Qu, Yixin Zhu, Andre Pradhana, Chenfanfu Jiang


** Build Instructions:

Step 1: Download and unzip mls-mpm88.zip (Link: http://bit.ly/mls-mpm88)
		Now you should have "mls-mpm88.cpp" and "taichi.h".

Step 2: Compile and run

* Linux:
	g++ mls-mpm88-explained.cpp -std=c++14 -g -lX11 -lpthread -O3 -o mls-mpm
	./mls-mpm


* Windows (MinGW):
	g++ mls-mpm88-explained.cpp -std=c++14 -lgdi32 -lpthread -O3 -o mls-mpm
	.\mls-mpm.exe


* Windows (Visual Studio 2017+):
  - Create an "Empty Project"
  - Use taichi.h as the only header, and mls-mpm88-explained.cpp as the only source
  - Change configuration to "Release" and "x64"
  - Press F5 to compile and run


* OS X:
	g++ mls-mpm88-explained.cpp -std=c++14 -framework Cocoa -lpthread -O3 -o mls-mpm
	./mls-mpm


** FAQ:
Q1: What does "1e-4_f" mean?
A1: The same as 1e-4f, a float precision real number.

Q2: What is "real"?
A2: real = float in this file.

Q3: What are the hex numbers like 0xED553B?
A3: They are RGB color values.
	The color scheme is borrowed from
	https://color.adobe.com/Copy-of-Copy-of-Core-color-theme-11449181/

Q4: How can I get higher-quality?
A4: Change n to 320; Change dt to 1e-5; Change E to 2e4;
	Change particle per cube from 500 to 8000 (Ln 72).
	After the change the whole animation takes ~3 minutes on my computer.

Q5: How to record the animation?
A5: Uncomment Ln 2 and 85 and create a folder named "tmp".
	The frames will be saved to "tmp/XXXXX.png".

	To get a video, you can use ffmpeg. If you already have taichi installed,
	you can simply go to the "tmp" folder and execute

	  ti video 60

	where 60 stands for 60 FPS. A file named "video.mp4" is what you want.

Q6: How is taichi.h generated?
A6: Please check out my #include <taichi> talk:
	http://taichi.graphics/wp-content/uploads/2018/11/include_taichi.pdf
	and the generation script:
	https://github.com/yuanming-hu/taichi/blob/master/misc/amalgamate.py
	You can regenerate it using `ti amal`, if you have taichi installed.

Questions go to yuanming _at_ mit.edu
							or https://github.com/yuanming-hu/taichi_mpm/issues.

													  Last Update: March 6, 2019
													  Version 1.5

----------------------------------------------------------------------------- */
