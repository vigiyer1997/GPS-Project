// Vignesh Iyer
// Alexandrea Jee
// Nathan Mower
// Receiver Program

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;


public class Receiver {
	// Index counter initially at 0.
	static int Index = 0; 
	static String[] input;
	static double[][] data, J;
	static double[] x0, B, xV, Sk;
	static double TV;
	static int M;
	static double c = Double.parseDouble("2.997924580000000000E+08");
	static double R = Double.parseDouble("6.367444500000000000E+06");
	static double S = Double.parseDouble("8.616408999999999651E+04");
	static ArrayList<String> values;
	static ArrayList<String> output;
	static NumberFormat formatter;
		
	// Read-Write methods. 
	//Read() method allows data to be interpreted and then executed in a Write() method.
	//Refer to CS 2420 assignments/notes for Read() Write() methods.
	static void ReadData() throws IOException {
		BufferedReader In = new BufferedReader(new InputStreamReader(System.in));
		String satelitte_lines;
		while((satelitte_lines = In.readLine()) != null) {
			for(String st : satelitte_lines.split(" ")) {
				values.add(st);
			}
		}
			int cap = values.size();
			input = new String[cap];
			for(int i = 0; i < cap; i++) { 
				input[i] = values.get(i); 
			}
			In.close();	
	}
	
	static void WriteData() {
		try(Writer writer = new BufferedWriter(new OutputStreamWriter
				(new FileOutputStream("receiver.log"),"UTF-8")))
		{
			for(int i = 0; i<output.size(); i++) { writer.write(output.get(i)+"\n"); }
			writer.close();
		}
		catch(IOException e) { 
			System.out.println("Throw exception. Halt Program."); 
			}
	}
	
	private static String[] RTD(double a) {
		//Converts Radians To Degrees, then converts to DMS (Refer to Homework 1 exercises)
		double 	b = a*180/Math.PI;
		if(a<0){b = -1*b;}
		int 	d = ((int) Math.floor(b));
		int 	m = ((int) Math.floor(60*(b-d)));
		double 	s = 60*(60*(b-d)-m);
		if(a<0)	{ return new String[]{""+d,""+m,""+s, ""+(-1)};	}
		return new String[]{""+d,""+m,""+s, ""+1};
	}
	
	static double TimePos(double[] xV) {
		//Recommended using m-1 equations of the form (13) and apply least squares approach. 
		//Refer to Exercise 10
		double values = data[0][0] + (1/c) * Math.sqrt(Math.pow(xV[0]-data[0][1], 2)+Math.pow(xV[1]-data[0][2], 2)+Math.pow(xV[2]-data[0][3], 2));
		return values;
	}
	
	// Difference between list indices a[i] and b[i]
	static double[] diff(double[] a, double[] b) {
		if(a.length==b.length)
		{
			double[] values = new double[a.length];
			for(int i = 0; i<a.length; i++) { 
				values[i]=a[i]-b[i]; 
			}
			return values;
		}
		return null;
	}
	
	static double norm(double[] a) {
		double sum = 0.0d;
		for(Double d : a){ sum = sum + d*d; }
		return Math.sqrt(sum);
	}
		
	static double[] Gradient_f(double[] x) {
		double[][] diffs = new double[M][3];
		for(int j = 0; j<M; j++) { 
			diffs[j] = diff(new double[]{data[j][1],data[j][2],data[j][3]},x); 
		}
		
		double[] N = new double[M];
		for(int j = 0; j<M; j++) { N[j] = norm(diffs[j]); }
		
		double[] A = new double[M-1];
		for(int j = 0; j<M-1; j++) { 
			A[j] = N[j+1]-N[j]-c*(data[j][0]-data[j+1][0]); 
		}
		
		double[][] XYZ = new double[3][M-1];
		for(int i = 0; i<3; i++)
		{
			for(int j = 0; j<M-1; j++)
			{
				XYZ[i][j] = diffs[j][i]/N[j] - diffs[j+1][i]/N[j+1];
			}
		}
		double[] values = new double[3];
		for(int j = 0; j<3; j++)
		{
			values[j] = 0.0d;
			for(int i = 0; i<M-1; i++) { values[j] += A[i]*XYZ[j][i]; }
			values[j] = 2.0d*values[j];
		}
		return values;
		
	}
	
	static double GradientF(int i, int j, double[] x)
	{
		// Solve by Newton's Method GradF=0
		// i-j returned element after Gradient F
		// Note: Positive Definite
		double values = 0.5d * (2.0d * x[j] - 2.0d * data[i][j+1]) / (Math.sqrt(Math.pow((x[0]-data[i][1]),2) + Math.pow((x[1]-data[i][2]),2) + Math.pow((x[2]-data[i][3]),2)));
		values -= 0.5d * (2.0d * x[j] - 2.0d * data[i+1][j+1]) / (Math.sqrt(Math.pow((x[0]-data[i+1][1]),2) + Math.pow((x[1]-data[i+1][2]),2) + Math.pow((x[2]-data[i+1][3]),2)));
		return values;
	}
	
	static double J_ij(int i, int j, double[] x)
	{
		//Jacobian Matrix (i-j indices) components
		double sum = 0.0d;
		for(int k = 0; k<M-1; k++) { 
			sum += GradientF(k, i, x) * GradientF(k, j, x); 
		}
		return 2.0d * sum;
	}
		
	static double[][] J_Matrix(double[] x)
	{
		//Jacobian Matrix 
		double[][] values = new double[3][3];
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				values[i][j] = J_ij(i, j, x);
			}
		}
		return values;
	}
	
	static double[] Matrix(double[][] A, double[] b) {
		//i-j breakdown
		double x = (A[0][1]*A[1][2]*b[2] - A[0][1]*A[2][2]*b[1] - A[0][2]*A[1][1]*b[2] + A[0][2]*A[2][1]*b[1] + A[1][1]*A[2][2]*b[0] - A[1][2]*A[2][1]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0]);
		double y = -1 * (A[0][0]*A[1][2]*b[2] - A[0][0]*A[2][2]*b[1] - A[0][2]*A[1][0]*b[2] + A[0][2]*A[2][0]*b[1] + A[1][0]*A[2][2]*b[0] - A[1][2]*A[2][0]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0]);
		double z = (A[0][0]*A[1][1]*b[2] - A[0][0]*A[2][1]*b[1] - A[0][1]*A[1][0]*b[2] + A[0][1]*A[2][0]*b[1] + A[1][0]*A[2][1]*b[0] - A[1][1]*A[2][0]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0]);
		return new double[]{x,y,z};
	}
	
	static double[] NewtonsMethod(double[] x, int d) {
		// Note: Have to calculate the gradient of the Hessian to solve GradF=0. Point xV where Grad F = 0.
		// Note: Partial derivatives of the gradient
		Sk = Matrix(J_Matrix(x), Gradient_f(x));
		double[] values = x.clone();
		for(int i = 0; i<3; i++) { values[i] -= Sk[i]; }
		if(Math.abs(diff(values,x)[0])<0.0000001d && Math.abs(diff(values,x)[1])<0.0000001d && Math.abs(diff(values,x)[2])<0.0000001d) { 
			return values; 
		}
		else if(d>9) { 
			return null; 
		}
		else
		{
			d++;
			return NewtonsMethod(values, d);
		}
	}
	
	static int Get_Bounds() {
		int m = 1;
		while(Index+1+5*m < input.length && Math.abs(Double.parseDouble(input[Index+1+5*(m)])-Double.parseDouble(input[Index+1+5*(m-1)]))<0.5d) { m++; }
		return m;
	}
	
	static double[][] ReadSatData() {
		int bounds = Get_Bounds();
		double[][] values = new double[bounds][4];
		// first input > last input
		for(int j = 0; j < bounds; j++) {
			// 5 values incremented through Index+4
			// Input 1 is vehicle Input 2 is broadcast
			values[j][0] = Double.parseDouble(input[Index + 1]);
			// [3,4,5] - (x,y,z) of the satellite broadcast timing
			values[j][1] = Double.parseDouble(input[Index + 2]);
			values[j][2] = Double.parseDouble(input[Index + 3]);
			values[j][3] = Double.parseDouble(input[Index + 4]);


			// Increment Index +5
			Index = Index + 5;
		}
		M = bounds;
		return values;
	}
	
	private static double[] R3(double alpha, double[] x) {
		return new double[] { Math.cos(alpha)*x[0] - Math.sin(alpha)*x[1], Math.sin(alpha)*x[0] + Math.cos(alpha)*x[1], x[2] };
	}
	
	private static String[] GeographicCoordinates(String[] a) {
		double t = Double.parseDouble(a[0]);
		double x1 = Double.parseDouble(a[1]);
		double y1 = Double.parseDouble(a[2]);
		double z1 = Double.parseDouble(a[3]);
		
		double[] xyz = R3(-2*Math.PI*t/S, new double[]{x1,y1,z1});
		double x = xyz[0];
		double y = xyz[1];
		double z = xyz[2];
		
		double psi;
		if(x*x+y*y == 0) {
			if(z>=0) { 
				psi = Math.PI/2;	
			}
			else { 
				psi = -1*Math.PI/2;	
			}
		}
		else { 
			psi = Math.atan2(z,Math.sqrt(x*x+y*y)); 
		}
		double lambda;
		if(x>0 && y>0) { 
			lambda = Math.atan2(y,x); 
		}
		else if(x < 0) { 
			lambda = Math.PI + Math.atan2(y,x); 
		}
		else { 
			lambda = 2*Math.PI + Math.atan2(y,x); 
		}
		lambda-=Math.PI;
		
		String[] Psi = RTD(psi);
		String[] Lambda = RTD(lambda);
		
		double h = Math.sqrt(x*x + y*y + z*z) - R;
		
		return new String[]{""+t, Psi[0], Psi[1], Psi[2], Psi[3], Lambda[0], Lambda[1], Lambda[2], Lambda[3], ""+h};
	}
		
	static void GetLocation() {
		// Calculates the location
		xV = NewtonsMethod(x0, 0);
		TV = TimePos(xV);
		String[] values = GeographicCoordinates(new String[]{""+TV, ""+xV[0], ""+xV[1], ""+xV[2]});
		String ret = "";
		for(int i = 0; i<10; i++)
		{
			ret = ret + formatter.format(Double.parseDouble(values[i])) + " ";
		}
		System.out.println(ret);
		output.add(ret);

	}
	
	public static void main(String[] args) throws IOException {
		//Main function; arbritary start point
		formatter = new DecimalFormat("#0.00");
		output = new ArrayList<String>();
		values = new ArrayList<String>();
		ReadData();
		x0 = new double[]{ 0, 0, 0 };
		while(Index<input.length)
		{
			data = ReadSatData(); // data is an m by n matrix
			GetLocation();
		}
		WriteData();
		
	}

}
