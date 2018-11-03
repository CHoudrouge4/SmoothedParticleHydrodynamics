import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner;

public class MeshConvertor {
	
	public static void convertToOFFfromMesh(String fileName) throws FileNotFoundException {
		Scanner scan = new Scanner(new File(fileName));
		PrintWriter print = new PrintWriter("out.off");
		int vert = 0, elem = 0;
		scan.nextLine();
		String str = "";
		String coord = ""; 
		while(scan.hasNextLine()) {
			String s = scan.nextLine();
			if(s.equals("elements"))
			{
				int size = scan.nextInt();
				for(int i = 0; i < size; ++i) {
					scan.nextInt(); 
					int type = scan.nextInt();
					if(type == 2) {
						str += "3 " + scan.nextInt() + " " + scan.nextInt() + " " + scan.nextInt() + '\n'; 		
						elem++;
					} else {
						scan.nextLine();
					}
				} 
			}
			if(s.equals("vertices")) {
				vert = scan.nextInt();
				while(!scan.hasNextDouble()) scan.nextLine();
				for(int i = 0; i < vert; ++i) {
					coord += scan.nextFloat() + " " + scan.nextFloat() + " " + scan.nextFloat() + '\n';
				}
			}
		}
		
		str = "OFF\n" +  vert + " " + elem + " " + 0 + "\n" + coord + '\n' + str;
		print.print(str);
		
		print.close();
		scan.close();
	}
}
