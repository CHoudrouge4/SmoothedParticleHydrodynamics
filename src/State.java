/*
 * Hussein Houdrouge
 */
public class State {
	
	int n;
	double m;
	// special coordinates
	double x[];
	double y[];
	double z[];
	
	// velocities values 
	double vx[];
	double vy[];
	double vz[];
	
	// half velocity values
	double vhx[];
	double vhy[];
	double vhz[];
	
	//acceleration
	double ax[];
	double ay[];
	double az[];
	
	double rho[];
	
	int height;
	int width;
	/*
	 *  n is the number of particles
	 *  height and width are the dimensions of the domain 
	 */
	public State(int n, int height, int width) {
		this.n = n;
		this.height = height;
		this.width = width;
		m   = 1;
		x   = new double[n];
		y   = new double[n];
		z   = new double[n];
		
		vx  = new double[n];
		vy  = new double[n];
		vz  = new double[n];
		
		vhx = new double[n];
		vhy = new double[n];
	    vhz = new double[n];
	    
		ax  = new double[n];
		ay  = new double[n];
		az  = new double[n];
		
		rho = new double[n];
		
		for(int i = 0; i < n; i++) {
			vx[i] = (Math.random() * 2 - 1)/1000;
		}
	}
}
