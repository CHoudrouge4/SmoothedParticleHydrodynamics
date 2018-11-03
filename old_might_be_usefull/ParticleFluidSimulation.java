import java.util.ArrayList;
import java.util.LinkedList;
public class ParticleFluidSimulation {
	int number_of_particle;
	int x_boundary;
	int left_boundary;
	int right_boundary;
	int up_boundary;
	int down_boundary;
	boolean ok = true;
	
	double p_m = 0.0002645833333333;
	double g = 9.8;
	
	LinkedList<Particle> lp;
	ArrayList<Integer>[][] grid;
	int h;
	int k = 30;
	int p0 = 20;
	double delta_t = 0.1;
	int hh;
	double poly = 1.27323954;//1.56668;	
	int raduis = 5;
	int h2, h4, h8;
	double C0, C1, C2;
	double mass = 1;
	double Cp;
	double rho0 = 1000;
	double rho02 = 2 * rho0;
	int mu = 3;
	double Cv = -40 * mu;
	
	public ParticleFluidSimulation(int n, int m, int step, int lb, int rb, int up, int db) {
		h = step;
		hh = 2 * 15;
		h2 = hh * hh;
		h4 = h2 * h2;
		h8 = h4 * h4;
		poly *= (1/h8);
		
		C0 = mass / (Math.PI * h4);
	    C1 = 4 * mass / (Math.PI * h2);
	    C2 = 4 * mass / (Math.PI * h8);
		Cp = 15 * k;
		
		number_of_particle = n * m;
		left_boundary = lb;
		right_boundary = rb;
		up_boundary = up;
		down_boundary = db;
		lp = new LinkedList<Particle>();
		int [] rgb = {100, 100, 100};
		int id = 0;
		for(int i = 0; i < n; ++i) {
			for(int j = 0; j < m; ++j) {
				double vx = (Math.random() * 2 - 1);
				Particle p = new Particle(lb + step * (i + 1), up + step * (j + 1) , 0, vx, 0, 0, 1, 1, 1,  g, 0, 0);
				p.id = id;
				id++;
				p.d = 1;
				lp.add(p);
			}
		}
		build_LookUp_Grid();
	}
	
	private int get_x_cell(double x) {
		double xx = x - left_boundary;
		return Math.max((int) Math.floor(xx / h) - 1, 0);
	}
	
	private int get_y_cell(double y) {
		double yy = y - up_boundary;
		return Math.max((int) Math.floor(yy / h) - 1, 0);
	}
	
	private void build_LookUp_Grid() {
		int y_size = Math.abs(up_boundary - down_boundary) / h;
		int x_size = Math.abs(left_boundary - right_boundary)/h;
		grid = new ArrayList[x_size][y_size];
		grid[0][0] = new ArrayList<Integer>(); // add another ArrayList object to [0,0]
		for(int i = 0; i < x_size; i++) {
			for(int j = 0; j < y_size; j++) {
				grid[i][j] = new ArrayList<Integer>();
			}
		}
		
		for(int i = 0; i < lp.size(); i++) {
			int index_x = get_x_cell(lp.get(i).x);
			int index_y = get_y_cell(lp.get(i).y);
			if(is_valid_ix(index_x) && is_valid_ix(index_y))
				grid[index_y][index_x].add(lp.get(i).id);
		}
	}
	
	public void compute_a_pressure(int i) {
		Particle p = lp.get(i);
		double rhoi = lp.get(i).d;
		p.ax = 0;
		p.ay = g;
		for(int j = i; j >= 0; j--) {
			double dx = p.x - lp.get(j).x;
			double dy = p.y - lp.get(j).y;
			double r2 = dx * dx + dy * dy;
			
			if(r2 < h2) {
				double rhoj = lp.get(j).d;
				
				double q = Math.sqrt(r2)/hh;
				System.out.println("q: " + q);
				
				double u = 1 - q;
				double w0 = C0 * u / (rhoi * rhoj);
				System.out.println("rhos " + rhoi * rhoj);
				
				double wp = w0 * Cp * (rhoi + rhoj - rho02) * u / q;
			    double wv = w0 * Cv;
			
			    System.out.println("w0 " + w0 + " " + "wp " + wp + " " +  "wv " +  wv); 
			    double dvx = p.vx - lp.get(j).vx;
			    double dvy = p.vy - lp.get(j).vy;
			    
			    p.ax += wp * dx + wv * dvx;
                p.ay += wp * dy + wv * dvy;
                lp.get(j).ax -= wp * dx + wv * dvx;
                lp.get(j).ay -= wp * dy + wv * dvy;
                p.ay *= -1;
                
			}
		}
	}

	public void update_x() {
		int n = lp.size();
		for(int i = 0; i < n; ++i) {
			Particle p = lp.get(i);
			p.vx = p.vx + p.ax * delta_t;
			p.vy = p.vy + p.ay * delta_t;
			System.out.println(p.ax);
			if(p.y + p.vy * delta_t > down_boundary || p.y + p.vy * delta_t < up_boundary ) {
				p.vy = -0.5* p.vy;	
				if(Math.abs(p.vy) <= 5) p.vy = 0;
			} else {
				p.y = p.y + p.vy * delta_t  +  0.5 * p.ay * delta_t * delta_t;
			}
			
			if(p.x + p.vx * delta_t < left_boundary || p.x + p.vx * delta_t > right_boundary ) {
				p.vx = -0.5*p.vx;
				if(Math.abs(p.vx) <= 5) p.vx = 0;
			} else {
				p.x  = p.x +  p.vx * delta_t +  0.5 * p.ax * delta_t * delta_t;
				p.ax = 0;
			}
		}
	}
	
	public void simulate(double t) {
		if(ok) {		
			compute_density(1);
			for(int i = 0; i < lp.size(); ++i) {
				lp.get(i).ax = 0;
				lp.get(i).ay = -g;
			}
			
			for(int i = lp.size() - 1; i >= 0; i--) {
				compute_a_pressure(i);
			} 
			
			update_x();
		}
		ok = false;
	}
	
	public double Delta_W(double d) {
		if(d > hh) return 0;
		return - 1 * 45/(Math.PI * h4) * (hh - d) * (hh - d);
	}
	
	public double W(double d) {
		if(d >= hh) return 0;
		System.out.println("doing computation");
		double term = (h2  -  d * d);
		return poly * Math.pow(term, 3);	
	}
	
	public double distance(int p, int q) {
		double px = lp.get(p).x;
		double py = lp.get(p).y;
		double qx = lp.get(q).x;
		double qy = lp.get(q).y;
		return (double) Math.sqrt((px - qx)*(px - qx) + (py - qy) * (py - qy));
	}
	
	public void compute_density(double mass) {
		double dx, dy, r2, z, rho_ij;
		double C1 = 4 * mass / (Math.PI * h2);
		double C2 = 4 * mass / (Math.PI * h8);
		
		for(int i = lp.size() - 1; i >= 0; i--)	lp.get(i).d = C1;
		
		for(int i = 0; i < lp.size(); ++i) {
			for(int j = i + 1; j < lp.size(); ++j) {
				dx = lp.get(i).x - lp.get(j).x;
				dy = lp.get(i).y - lp.get(j).y;
				r2 = dx * dx + dy * dy;
				z = h2 - r2;
			
				if(z > 0) {
					rho_ij = C2 * z * z * z;
					lp.get(i).d += rho_ij;
					lp.get(j).d += rho_ij;
				}
			}
		}
	}
	
	public void pressure(int p) {
		lp.get(p).p = k * (lp.get(p).d - p0);
	}
	
	private boolean is_valid_ix(int index_x) {
		return (index_x >= 0 && index_x < grid[0].length);
	}
	
	private boolean is_valid_iy(int index_y) {
		return (index_y >= 0 && index_y < grid.length);
	}
	
	public LinkedList<Integer> neighbors(int p) {
		LinkedList<Integer> list = new LinkedList<>();
		int index_x = get_x_cell(lp.get(p).x);
		int index_y = get_y_cell(lp.get(p).y);
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
	
	public void simulate_one() {
		for(int i = 0; i < lp.size(); ++i) {
			Particle p = lp.get(i);
			p.ax = 0;
			p.ay = 9;
			
			p.vx = p.vx + p.ax * delta_t;
			p.vy = p.vy + p.ay * delta_t;
	
			if(p.y + p.vy * delta_t > down_boundary) {
				p.vy = -0.5* p.vy;	
				if(Math.abs(p.vy) <= 5) p.vy = 0;
			}
			
			p.x = p.x + p.vx * delta_t;
			p.y = p.y + p.vy * delta_t + 0.5 * p.ax * delta_t *delta_t;
		}
	}
}
