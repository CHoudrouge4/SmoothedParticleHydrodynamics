import processing.core.*;

public class Viewer extends PApplet {
	
	ParticleFluidSimulation pfs;
	State s;
	ArcBall arcball;
	public void setup() {		
		  size(600, 600,P2D);
		  arcball = new ArcBall(this);
		   
		  ellipseMode(RADIUS);
		  s = new State(1000, HEIGHT, WIDTH);
		  pfs = new ParticleFluidSimulation(s);
		 
		  pfs.init();
		
		 // LinkedList<Integer> ll = pfs.neighbors(0);
	}
	
	public void draw() {
		  background(0);
		  this.lights();
		  directionalLight(101, 204, 255, -1,  0, 0);
		  directionalLight(51,  102, 126,  0, -1, 0);
		  directionalLight(51,  102, 126,  0,  0, -1);
		  directionalLight(102,  50, 126,  1,  0, 0);
		  directionalLight(51,   50, 102,  0,  1, 0);
		  directionalLight(51,   50, 102,  0,  0, 1);  
		  
		// translate(0, 0, 0);
			  
		 /*  this.camera(0.70f * width, 0.35f * width, 1.20f * width, width/2, width/2, 0.0f, 0.0f, 1.0f, 0.0f);	
			 translate(50, 50, 0);
			 rotateX(-3 *PI/3);
			 rotateY(3 * PI/3);
			 rotateZ(PI/2);
		  */
		
		 this.strokeWeight(1);
				
		 for(int i = 0; i < s.n; i++) {
			 this.ellipse((float)(s.x[i]) * height, (float)(s.y[i]) * height, 3f , 3f);
			 //noStroke();
			 //translate((float)(s.x[i]) * height, (float)(s.y[i]) * height, (float) (s.z[i]) * height);
			 //sphere(2.5f);
			 //translate(-1 * (float)(s.x[i]) * height, -1 * (float)(s.y[i]) * height, (float) (s.z[i]) * height);
		 }
		
		 //this.point(width/2, height/2, 0);
		 
		 pfs.simulate();
	}
		
	public void keyPressed(){
		if(key == 'i') zoomIn();
		if(key == 'o') zoomOut();
	}
	
	  public void zoomIn() {
		  System.out.println("zoom in");
		  this.arcball.scaleFactor=this.arcball.scaleFactor*1.1f;
	  }

	  public void zoomOut() {
		  System.out.println("zoom out");
		  this.arcball.scaleFactor= this.arcball.scaleFactor*0.9f;
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