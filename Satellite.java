// Vignesh Iyer
// Alexandrea Jee
// Nathan Mower
// Satellite Program

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Scanner;



public class Satellite {
	static double Pi, c, R, S;
	static double [][] V;
	static double [] Xv, Xs;
	static double tv, ts;
	static String path;
	static ArrayList <String> output;
	
	
	private static double DTR(String[] a) {
		return 2*Pi*(Integer.parseInt(a[0])/360.0d + Integer.parseInt(a[1])/(360.0d*60) + Double.parseDouble(a[2])/(360*60*60));
	}
	
	private static double[] CartesianCoordinates(String[] a) {
		double 	t = Double.parseDouble(a[0]);
		double 	theta = Integer.parseInt(a[4])*DTR(new String[]{a[1],a[2],a[3]});
		double 	phi	= Integer.parseInt(a[8])*DTR(new String[]{a[5],a[6],a[7]});
		double 	h = Double.parseDouble(a[9]);
		double 	x = (R + h) * Math.cos(theta) * Math.cos(phi);
		double 	y = (R + h) * Math.cos(theta) * Math.sin(phi);
		double 	z = (R + h) * Math.sin(theta);
		double  alpha = (2*Pi*t)/S;
		
		return new double[]{Math.cos(alpha)*x-Math.sin(alpha)*y,Math.sin(alpha)*x+Math.cos(alpha)*y,z,t};
	}
	
	private static double f(int v, double t) {
		double value = -1*c*c*(tv-t)*(tv-t);
		for(int i = 0; i<3; i++)
		{
			value += (Xs(v,t)[i]-Xv[i])*(Xs(v,t)[i]-Xv[i]);
		}
		return value;
	}
		
	private static double f1(int v, double t) {
		double value = 2*c*c*(tv-t);
		for(int i = 0; i<3; i++)
		{
			value += 2*(Xs(v,t)[i]-Xv[i])*Xs1(v,t)[i];
		}
		return value;
	}
	
	private static double NewtonsMethod(int v, double Tk, int depth) {
		double value = Tk;
		value = Tk-(f(v,Tk)/f1(v,Tk));
		if(value-Tk < 0.01/c) { 
			return value; 
		}
		else if(depth>=9) { 
			return -1; 
		}
		else { 
			return NewtonsMethod(v, value, depth+1); 
		}
	}
	
	private static double Get_ts(int v) {
		double t0 = tv;
		double norm = 0.0d;
		for(int i = 0; i < 3; i++)
		{
			norm += Math.pow(Xs(v,tv)[i],2);
		}
		t0 = t0 - Math.sqrt(norm)/c;
		
		return NewtonsMethod(v,t0,0);
	}
	
	private static double[] Xs(int v, double t) {
		double x = (R+V[v][7])*(V[v][0]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][3]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		double y = (R+V[v][7])*(V[v][1]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][4]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		double z = (R+V[v][7])*(V[v][2]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][5]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		return new double[]{x,y,z,t};
	}
	
	private static double[] Xs1(int v, double t) {
		double x = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][3]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][0]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		double y = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][4]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][1]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		double z = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][5]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][2]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		return new double[]{x,y,z,t};
	}
	
	// Compare xTs>xTx (Check notes)
	private static boolean[] above_Horizon(double[] a) {
		boolean[] value = new boolean[24];
		double[] satXYZ;
		for(int i = 0; i < 24; i++)
		{
			satXYZ = Xs(i, a[3]);
			value[i] = 2*a[0]*(satXYZ[0]-a[0]) + 2*a[1]*(satXYZ[1]-a[1]) + 2*a[2]*(satXYZ[2]-a[2]) > 0;
		}
		return value;
	}
	
	private static void Initialize()throws IOException {
		V = new double[24][9];
		output = new ArrayList<String>();
		String temp = "";
		Scanner scan = new Scanner(new File("data.dat")).useDelimiter("\n"); //disregard closable value, efficiency not counted.
		scan.useDelimiter("\n");
		String[] l = new String[220];
		int counter = 0;
		while(scan.hasNext())
		{
			temp = scan.next();
			l[counter] = temp;
			counter++;
		}
		scan.close();
		Pi = Double.parseDouble(l[0].substring(0, 27));
		c = Double.parseDouble(l[1].substring(0, 27));
		R = Double.parseDouble(l[2].substring(0, 27));
		S = Double.parseDouble(l[3].substring(0, 27));
		
		counter = 4;
		for(int v = 0; v < 24; v++)
		{
			for(int i = 0; i < 9; i++)
			{
				V[v][i] = Double.parseDouble(l[counter].substring(0, 27));
				counter++;
			}
		}
	}
	
	private static void ProcessLine(String s) {
		String[] split = s.split(" ");
		double[] xyz = CartesianCoordinates(split);
		tv = Double.parseDouble(split[0]);
		Xv = xyz;
		boolean[] B = above_Horizon(xyz);
		for(int i = 0; i < B.length; i++)
		{
			if(B[i])
			{
				ts = Get_ts(i);
				Xs = Xs(i, ts);
				String ret = i+" "+ts+" "+Xs[0]+" "+Xs[1]+" "+Xs[2];
				System.out.println(ret);
				output.add(ret);
				if(ts>200 && ts<203)
				{
					System.out.println(":   ");
					output.add(":   ");
				}
			}
		}
	}

	private static void ReadDataFile() throws IOException {
		BufferedReader In = new BufferedReader(new InputStreamReader(System.in));
		String s;
		while((s = In.readLine()) != null) { 
			ProcessLine(s); 
		}
	}
	
	static void WriteDataFile() {
		try(Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("satellite.log"),"UTF-8")))
		{
			for(int i = 0; i<output.size(); i++) { 
				writer.write(output.get(i)+"\n"); 
			}
			writer.close();
		}
		catch(IOException e) { 
			System.out.println("Catch exception and halt program"); 
		}
	}
	
	public static void main(String[] args) throws IOException {
		path = ".\\";
		
		ts = 0.0d;
		Xs = null;
		Initialize();
		ReadDataFile();
		WriteDataFile();
	}

	
	

}
