
public class Particle {
	int id;
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double M;
	double d; // density
	double p; // pressure
	int [] rgb;
	double fx;
	double fy;
	double fz;
	double ax;
	double ay;
	double az;
	
	public Particle(double x, double y, double z, double vx, double vy, double vz, double M, double d, double p, double fx, double fy, double fz) {
		this.x = x;
		this.y = y;
		this.z = z;
		this.vx = vx;
		this.vy = vy;
		this.vz = vz;
		this.M = M;
		this.d = d;
		this.p = p;
		this.fx = fx;
		this.fy = fy;
		this.fz = fz;
	}
}
