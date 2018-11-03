import java.util.ArrayList;
import java.util.LinkedList;

/*
 * Hussein Houdrouge 
 * - 
 * This code is based on the following instructions 
 * http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf
 */

public class ParticleFluidSimulation {

	double h = 0.02;   				// Raduis of the particle
	double h2 = h * h;
	double h4 = h2 * h2;
	double h8 = h4 * h4;
	double dt = 10e-4;		   // delta t 
	double rho0 = 800;		  // Initial Concentration 
	double k = 30;           // Bulk modulus
	double mu = 2;			// Viscosity
	double g = -10;        // Gravity 
	State s;
	double pi = Math.PI;// 3.141526f;
	double cp = 15 * k;
	double cv = -40 * mu;
	double c0, c1, c2;
	ArrayList<Integer>[][] grid;
	
	
	public ParticleFluidSimulation(State s) {
		this.s = s;
	}
	
	/*
	 *  This function compute the density 
	 *  of each particle
	 */
	private void compute_density() {
		int n = s.n;
		c1 = 4 * s.m / (pi * h2);
		c2 = 4 * s.m / (pi * h8);
		//add
		for(int i = 0; i < n; i++) s.rho[i] = c1;
		
		for(int i = 0; i < n; i++) {
			for(int j = i + 1; j < n; j++) {
				double dx = s.x[i] - s.x[j];
				double dy = s.y[i] - s.y[j];
				double r2 = dx * dx + dy * dy;
				double z = h2 - r2;
				if(z > 0) {
					double rho_ij = c2 * z * z * z;
					s.rho[i] += rho_ij;
					s.rho[j] += rho_ij;
				}
			}
		}
	}
	
	/*
	 *  This function compute the acceleration for each particle
	 * 
	 */
	public void compute_accel() {
		int n = s.n;
		// compute density
		compute_density();
		
		for(int i = 0; i < n; ++i) {
			s.ax[i] = 0;
			s.ay[i] = -g;
			s.az[i] = 0;
		}
		
		for(int i = 0; i < n; i++) {
			double rhoi = s.rho[i];
			for(int j = i + 1; j < n; j++) {
				double dx = s.x[i] - s.x[j];
				double dy = s.y[i] - s.y[j];
				double dz = s.z[i] - s.z[j];
				double r2 = dx * dx + dy * dy + dz * dz;
				
				if(r2 < h2) {
					double rhoj = s.rho[j];
					double q = Math.sqrt((double)r2)/h;
					double u = 1.0 - q;
					double w0 = c0 * u / (rhoi * rhoj);
					double wp = w0 * cp * (rhoi + rhoj - 2 * rho0) * u /q;
					double wv = w0 * cv;
	
					double dvx = s.vx[i] - s.vx[j];
					double dvy = s.vy[i] - s.vy[j];
					double dvz = s.vz[i] - s.vz[j];
					
					double valx = (wp * dx + wv * dvx);
					double valy = (wp * dy + wv * dvy);
					double valz = (wp * dz + wv * dvz); 
					s.ax[i] += valx;
					s.ay[i] += valy;
					s.az[i] += valz;
		
					s.ax[j] -= valx;
					s.ay[j] -= valy;
					s.az[j] -= valz;
				}
			}
		}
	}
	
	
	/*
	 * 
	 *  This function implement leapfrog at each iteration 
	 *  
	 */
	public void leapfrog_step() {
		int n = s.n;
		for(int i = 0; i < n; ++i) {
			s.vhx[i] += s.ax[i] * dt;
			s.vhy[i] += s.ay[i] * dt;
			s.vhz[i] += s.az[i] * dt;
			
			s.vx[i] = s.vhx[i] + s.ax[i] * dt/2;
			s.vy[i] = s.vhy[i] + s.ay[i] * dt/2;
			s.vz[i] = s.vhz[i] + s.az[i] * dt/2;
		
			s.x[i] += s.vhx[i] * dt;
			s.y[i] += s.vhy[i] * dt;
			s.z[i] += s.vhz[i] * dt;
		}
		
		boundary_reflection();
	}
	
	/*
	 * 
	 * This function implement the first step in leapfrog
	 * compute V0, and V1/2 for each particles
	 */
	public void leapfrog_start() {
		int n = s.n;
		for(int i = 0; i < n; i++) {
			s.vhx[i] = s.vx[i] + s.ax[i] * dt/2;
			s.vhy[i] = s.vy[i] + s.ay[i] * dt/2;
			s.vhz[i] = s.vz[i] + s.az[i] * dt/2;
			
			s.vx[i] = s.ax[i] * dt/2;
			s.vy[i] = s.ay[i] * dt/2;
			s.vz[i] = s.az[i] * dt/2;
		
			s.x[i] += s.vhx[i] * dt;
			s.y[i] += s.vhy[i] * dt;
			s.z[i] += s.vhz[i] * dt;
		}
		
		boundary_reflection();	
	}
	
	/*
	 * - reflect the particle by changing the speed.
	 *   (Multiplying by -a where a between 0 and 1)
	 * 
	 */
	private void reflection(int which, double barrier, int i) {
		if(which == 1) {
			s.vy[i]  = -0.5* s.vy[i];
			s.vhy[i] = -0.5*s.vhy[i];
			if(Math.abs(s.vy[i]) <= 0.05) s.vy[i] = 0;
			s.y[i] += s.vy[i] * dt;
		}

		if(which == 0) {
			s.vx[i]  = -0.5 * s.vx[i];
			s.vhx[i] = -0.5* s.vhx[i];
			if(Math.abs(s.vx[i]) <= 5) s.vx[i] = 0;
			s.x[i] += s.vx[i] *dt;
		} 
	
		if(which == 2) {
			s.vz[i] = -0.5 * s.vz[i];
			s.vhz[i] = -0.5 * s.vhz[i];
			if(Math.abs(s.vz[i]) <= 5) s.vz[i] = 0;
			s.z[i] += s.vz[i] *dt;
		}
		
	}
	
	private void boundary_reflection() {
		final double xmin = 0.0;
		final double xmax = 1.0;
		final double ymin = 0.0;
		final double ymax = 1.0;
		final double zmin = 0.0;
		final double zmax = 1.0;
		
		int n = s.n;
		for(int i = 0; i < n; ++i) {
			if(s.x[i] < xmin) reflection(0, xmin, i);
			if(s.x[i] > xmax) reflection(0, xmax, i);
			if(s.y[i] < ymin) reflection(1, ymin, i);
			if(s.y[i] > ymax) reflection(1, ymax, i);
			//if(s.z[i] < zmin) damp_reflect(2, zmin, i);
			//if(s.z[i] > zmax) damp_reflect(2, zmax, i);
		}
	}
	
	/*
	 * This function place the particle in a rectangular shape
	 * of length y2 - y1 
	 */
	private void place_particles_plane(double y1, double y2, double zp, int ni, int nf) {
		double r = h;
		double xp = h + 0.04;
		double yp = y1;
		for(int i = ni; i < nf; i++) {
			s.x[i] = xp;
			s.y[i] = yp;
			s.z[i] = zp;
			yp += r;
			if(yp > y2) {
					yp = y1;
					xp += r;	
			}
		}
	}
	
	public void place_particles(double y1, double y2, double z2) {
		double zp = h * 0.5 + 0.01;
		
		int nb_of_pp = s.n;
		
		int starts = 0;
		int end = nb_of_pp;
		//for(int i = 0 ; i < 200; i++) {
			place_particles_plane(y1, y2, zp, starts , end);
		//	starts += nb_of_pp;
		//	end += nb_of_pp;
			//zp += h;
		//}
	}
	
	private void normalize_mass() {
		s.m = 1;
		compute_density();
		//double rho0 = this.rho0;
		double rho2s = 0;
		double rhos = 0;
		for(int i = 0; i < s.n; i++) {
			rho2s += s.rho[i] * s.rho[i];
			rhos  += s.rho[i];
 		}
		s.m = (rho0 * rhos)/rho2s;
		c0 = s.m /(pi * h4);
		c1 = 4 * s.m / (pi * h2);
		c2 = 4 * s.m / (pi * h8);
	}
	
	/*
	 * 
	 * This function is one step in the simulation
	 */
	public void simulate() {
		compute_accel();
		leapfrog_step();
	}
	
	/*
	 * This function is the first step in the  
	 * simulation
	 */
	public void init() {
		place_particles(0.05, s.height / s.width - 0.01, s.height / s.width - 0.01);
		normalize_mass();
		compute_accel();
		leapfrog_start();
	}
	
	private boolean is_valid_ix(int index_x) {
		return (index_x >= 0 && index_x < grid[0].length);
	}
	
	private boolean is_valid_iy(int index_y) {
		return (index_y >= 0 && index_y < grid.length);
	}
	
	private int get_x_cell(double x) {
		double xx = x - 0;
		return Math.max((int) Math.floor(xx / h) - 1, 0);
	}
	
	private int get_y_cell(double y) {
		double yy = y - 0;
		return Math.max((int) Math.floor(yy / h) - 1, 0);
	}
	
	public LinkedList<Integer> neighbors(int p) {
		LinkedList<Integer> list = new LinkedList<>();
		int index_x = get_x_cell(s.x[p]);
		int index_y = get_y_cell(s.y[p]);
		if(is_valid_ix(index_x + 1) && is_valid_iy(index_y)) list.addAll(grid[index_y][index_x + 1]);
		if(is_valid_ix(index_x - 1) && is_valid_iy(index_y)) list.addAll(grid[index_y][index_x - 1]);
		if(is_valid_iy(index_y + 1) && is_valid_ix(index_x)) list.addAll(grid[index_y + 1][index_x]);
		if(is_valid_ix(index_y - 1) && is_valid_ix(index_x)) list.addAll(grid[index_y - 1][index_x]);
		if(is_valid_ix(index_x + 1) && is_valid_iy(index_y + 1)) list.addAll(grid[index_y + 1][index_x + 1]);
		if(is_valid_ix(index_x + 1) && is_valid_iy(index_y - 1)) list.addAll(grid[index_y - 1][index_x + 1]);
		if(is_valid_ix(index_x - 1) && is_valid_iy(index_y + 1)) list.addAll(grid[index_x + 1][index_x - 1]);
		if(is_valid_ix(index_x - 1) && is_valid_iy(index_y - 1)) list.addAll(grid[index_y - 1][index_x - 1]);
		if(is_valid_ix(index_x)     && is_valid_iy(index_y)) list.addAll(grid[index_y][index_x]);
		return list;
	}
	
	private void build_LookUp_Grid() {
		int y_size =(int) (1 / h);
		int x_size = (int) (1 / h);
		grid = new ArrayList[x_size][y_size];
		grid[0][0] = new ArrayList<Integer>(); // add another ArrayList object to [0,0]
		for(int i = 0; i < x_size; i++) {
			for(int j = 0; j < y_size; j++) {
				grid[i][j] = new ArrayList<Integer>();
			}
		}
		
		for(int i = 0; i < s.n; i++) {
			int index_x = get_x_cell(s.x[i]);
			int index_y = get_y_cell(s.y[i]);
			if(is_valid_ix(index_x) && is_valid_ix(index_y))
				grid[index_y][index_x].add(i);
		}
	}
	
}