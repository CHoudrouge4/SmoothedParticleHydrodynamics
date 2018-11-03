import processing.core.*;

public class Viewer extends PApplet {
	
	
	ParticleFluidSimulation pfs;
	public void setup() {		
		  size(1000, 600,P3D);
		  ArcBall arcball = new ArcBall(this); 
		  ellipseMode(RADIUS);
		  pfs = new ParticleFluidSimulation(10, 10, 15, 100, 800, 100, 500);
		 // LinkedList<Integer> ll = pfs.neighbors(0);
	}
	
	public void draw() {
		  background(0);
		  //this.lights();
		  directionalLight(101, 204, 255, -1,  0, 0);
		  directionalLight(51,  102, 126,  0, -1, 0);
		  directionalLight(51,  102, 126,  0,  0, -1);
		  directionalLight(102,  50, 126,  1,  0, 0);
		  directionalLight(51,   50, 102,  0,  1, 0);
		  directionalLight(51,   50, 102,  0,  0, 1);
		  
		//  this.rect(100, 100, 700, 500);
		  
		 // this.ellipse(100, 100, 10, 10);
		  pfs.simulate(0.001);
		  for(int i = 0; i < pfs.lp.size(); ++i) {
			  double x = pfs.lp.get(i).x;
			  double y = pfs.lp.get(i).y;
			  double z = pfs.lp.get(i).z;
			  if(x < 800)
				  this.ellipse((float ) x, (float) y, 5, 5);
		  } 
		  
		  
		 translate(width/2.f, height/2.f, -1 * height/2.f);
		 this.strokeWeight(1);
		// stroke(150, 150, 150); 
	}
		
	public void keyPressed(){
	
	}
		
	/**
	 * For running the PApplet as Java application
	 */
	public static void main(String args[]) {
		//PApplet pa=new MeshViewer();
		//pa.setSize(400, 400);
		PApplet.main(new String[] { "Viewer" });
	}
}