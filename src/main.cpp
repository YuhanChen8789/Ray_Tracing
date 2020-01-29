// C++ include
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;


double deta_ans(Vector3d ray_origin, Vector3d ray_direction, Vector3d center, double sphere_radius)
{
            double a = ray_direction.dot(ray_direction);
            double b = 2 * ray_direction. dot(ray_origin - center);
            double c = (ray_origin - center).dot(ray_origin - center) - sphere_radius * sphere_radius;
            double deta = pow(b, 2) - 4 * a * c;

            // The ray hit the sphere with one or two point, compute the exact intersection point
	    if(deta >= 0)
		{
              //  *deta_true = true;
		double t1 = ((-1)*b + sqrt(deta)) / 2 * a;
                double t = ((-1)*b - sqrt(deta)) / 2 * a;
                //if(abs(t1) < abs(t))
                //t = t1;
                return t;
		}
            else{
              // *deta_true = false;
                return 100.0;
	         }
            return 0.0;
}

double diffuse(Vector3d ray_origin, Vector3d ray_direction, Vector3d center, double t, Vector3d light_position)
{
                   Vector3d ray_intersection = ray_origin + t * ray_direction;

                   // Compute normal at the intersection point
                   Vector3d ray_normal = (ray_intersection - center).normalized();
                   
                   return (light_position-ray_intersection).normalized().dot(ray_normal);
                   }



void part1()
{
    std::cout << "Part 1: Simple ray tracer, three sphere with orthographic projection" << std::endl;

    const std::string filename("part1.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/C.cols(),0,0);
    Vector3d y_displacement(0,-2.0/C.rows(),0);

    //center of a circle
    Vector3d center1(0.6,0.6,1);
    Vector3d center2(-0.5,0.5,1);
    Vector3d center3(0,-0.5,1);
    
    // Single light source
    const Vector3d light_position(0,0,1);

    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            
            const double sphere_radius1 = 0.25;
            const double sphere_radius2 = 0.35;
            const double sphere_radius3 = 0.40;

	    // Prepare the ray
	    
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

            //Intersect with the sphere
           double t1 = deta_ans(ray_origin, ray_direction, center1, sphere_radius1);
           double t2 = deta_ans(ray_origin, ray_direction, center2, sphere_radius2);
           double t3 = deta_ans(ray_origin, ray_direction, center3, sphere_radius3);
           
	   //find the closest sphere
           double min_t = 100;
           int index = 0;
           if(t1 < min_t)
           {
		min_t = t1;
                index = 1;
           }
	   if(t2 < min_t)
           {
		min_t = t2;
                index = 2;
           }
	   if(t3 < min_t)
           {
		min_t = t3;
                index = 3;
           }
          // Simple diffuse model
          if(index == 1)
          {
	  C(i,j) = diffuse(ray_origin, ray_direction, center1, t1, light_position);
	  // Clamp to zero
          C(i,j) = max(C(i,j),0.);
          // Disable the alpha mask for this pixel
          A(i,j) = 1;}

          if(index == 2)
          {C(i,j) = diffuse(ray_origin, ray_direction, center2, t2, light_position);
	  // Clamp to zero
           C(i,j) = max(C(i,j),0.);
           // Disable the alpha mask for this pixel
           A(i,j) = 1;}
         
          if(index == 3)
          {C(i,j) = diffuse(ray_origin, ray_direction, center3, t3, light_position);
	  // Clamp to zero
          C(i,j) = max(C(i,j),0.);
          // Disable the alpha mask for this pixel
          A(i,j) = 1;}
                              
           
      }
           
    }

    // Save to png
    write_matrix_to_png(C,C,C,A,filename);

}

double specular(Vector3d ray_origin, Vector3d ray_direction, Vector3d center, double t, Vector3d light_position, Vector3d origin)
{

		Vector3d ray_intersection = ray_origin + t * ray_direction;
		// Compute normal at the intersection point
                Vector3d ray_normal = (ray_intersection - center).normalized();
		Vector3d v = (-1) * (ray_direction).normalized();
                Vector3d l = (light_position - ray_intersection).normalized();
                Vector3d h = (v+l).normalized();
                return ray_normal.dot(h);
} 
double _ambient = 0.1;

void part2()
{
    std::cout << "Part 2: Simple ray tracer, three sphere with orthographic projection" << std::endl;

    const std::string filename("part2.png");
    MatrixXd R = MatrixXd::Zero(800,800); 
    MatrixXd G = MatrixXd::Zero(800,800);
    MatrixXd B = MatrixXd::Zero(800,800);// Store the color
    MatrixXd A = MatrixXd::Zero(800,800);// Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/R.cols(),0,0);
    Vector3d y_displacement(0,-2.0/R.rows(),0);
    //Vector3d origin1(-0.5,0.5,1);

    //center of a circle
    Vector3d center1(0.6,0.6,1);
    Vector3d center2(-0.5,0.5,1);
    Vector3d center3(0,-0.4,1);
    
    // two light sources
    const Vector3d light_position(0.4, 0,1.4);
    const Vector3d light_position1(-1,-1,1.4);

    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {
            
            const double sphere_radius1 = 0.25;
            const double sphere_radius2 = 0.35;
            const double sphere_radius3 = 0.3;

	    // Prepare the ray
	    
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

            // Intersect with the sphere
          
           double t1 = deta_ans(ray_origin, ray_direction, center1, sphere_radius1);
           double t2 = deta_ans(ray_origin, ray_direction, center2, sphere_radius2);
           double t3 = deta_ans(ray_origin, ray_direction, center3, sphere_radius3);
	  
	   //find the closest sphere
           
           double min_t = 100;
           int index = 0;
           if(t1 < min_t )
           {
		min_t = t1;
                index = 1;
           }
	   if(t2 < min_t)
           {
		min_t = t2;
                index = 2;
           }
	   if(t3 < min_t)
           {
		min_t = t3;
                index = 3;
           }

          
         // Final Shading model
          if(index == 1)
          {
	  double _diffuse, _diffuse1;
	  _diffuse = diffuse(ray_origin, ray_direction, center1, t1, light_position);
          _diffuse = max(_diffuse, 0.); 
          _diffuse1 = diffuse(ray_origin, ray_direction, center1, t1, light_position1);
          _diffuse1 = max(_diffuse1, 0.); 
          _diffuse = _diffuse + _diffuse1;
          double _specular;
          _specular = specular(ray_origin, ray_direction, center1, t1, light_position, origin);
          _specular = pow(max(_specular, 0.),1000);
          R(i,j) =0*( _diffuse + _ambient + _specular);
          G(i,j) =0*( _diffuse + _ambient + _specular);
          B(i,j) =1*(_diffuse + _ambient);
          // Disable the alpha mask for this pixel
          A(i,j) = 1;


          }

          if(index == 2)
          {
          double _diffuse, _diffuse1;
	  _diffuse = diffuse(ray_origin, ray_direction, center2, t2, light_position);
          _diffuse = max(_diffuse, 0.); 
          _diffuse1 = diffuse(ray_origin, ray_direction, center2, t2, light_position1);
          _diffuse1 = max(_diffuse1, 0.); 
         _diffuse = _diffuse + _diffuse1;
          double _specular, _specular1;
          _specular = specular(ray_origin, ray_direction, center2, t2, light_position, origin);
          _specular = pow(max(_specular, 0.),200);
          _specular1 = specular(ray_origin, ray_direction, center2, t2, light_position1, origin);
          _specular1 = pow(max(_specular1, 0.),200);
          _specular = (_specular + _specular1);
         // cout << "  " << _specular;
          R(i,j) =0*( _diffuse + _ambient + _specular);
          G(i,j) =1*( _diffuse + _ambient );
         // cout <<  G(i,j) << " ";
          B(i,j) =0*(_diffuse + _ambient + _specular);

           // Disable the alpha mask for this pixel
           A(i,j) = 1;}
         
          if(index == 3)
          {
           double _diffuse, _diffuse1;
	  _diffuse = diffuse(ray_origin, ray_direction, center3, t3, light_position);
          _diffuse = max(_diffuse, 0.); 
          _diffuse1 = diffuse(ray_origin, ray_direction, center3, t3, light_position1);
          _diffuse1 = max(_diffuse1, 0.);
          _diffuse = _diffuse + _diffuse1; 
          double _specular, _specular1;
          //Vector3d origin1(-1,5,5);
          _specular = specular(ray_origin, ray_direction, center3, t3, light_position, ray_origin);
          _specular = pow(max(_specular, 0.),500);
          _specular1 = specular(ray_origin, ray_direction, center3, t3, light_position1, ray_origin);
          _specular1 = pow(max(_specular1, 0.),500);
          _specular =  (_specular + _specular1);

          R(i,j) =1*( _diffuse + _ambient+ _specular);
          G(i,j) =0*( _diffuse + _specular + _ambient);
          B(i,j) =0*(_diffuse + _specular + _ambient);


          // Disable the alpha mask for this pixel
          A(i,j) = 1;}
                              
      }
           
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);

}

void part3()
{
    std::cout << "Part 3: Simple ray tracer, one sphere with perspective projection" << std::endl;

    const std::string filename("part3.png");
    MatrixXd R = MatrixXd::Zero(800,800); 
    MatrixXd G = MatrixXd::Zero(800,800);
    MatrixXd B = MatrixXd::Zero(800,800);// Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is persepective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(0,0,1.85);
    Vector3d x_displacement(2.0/R.cols(),0,0);
    Vector3d y_displacement(0,-2.0/R.rows(),0);
    double distance = 0.86;
    Vector3d w(0,0,1);
    Vector3d w2(-1,1,0);

    //center of a circle
    Vector3d center1(0.6,0.6,1);
    Vector3d center2(-0.5,0.5,1);
    Vector3d center3(0,-0.4,1);
    
    // two light sources
    const Vector3d light_position(0, 0,1.4);
    const Vector3d light_position1(-1,-1,1.4);

    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {
            
            const double sphere_radius1 = 0.25;
            const double sphere_radius2 = 0.35;
            const double sphere_radius3 = 0.3;

	    // Prepare the ray
	    
            Vector3d ray_origin = origin;
            Vector3d ray_direction = (-1 * distance * w + w2 +double(i)*x_displacement + double(j)*y_displacement).transpose();

            // Intersect with the sphere
         
          
           double t1 = deta_ans(ray_origin, ray_direction, center1, sphere_radius1);
           double t2 = deta_ans(ray_origin, ray_direction, center2, sphere_radius2);
           double t3 = deta_ans(ray_origin, ray_direction, center3, sphere_radius3);

          //find the closest sphere
           double min_t = 100;
           int index;
           if(t1 < min_t)
           {
		min_t = t1;
                index = 1;
           }
	   if(t2 < min_t)
           {
		min_t = t2;
                index = 2;
           }
	   if(t3 < min_t)
           {
		min_t = t3;
                index = 3;
           }

         // Final Shading model
          if(index == 1)
          {
	  double _diffuse, _diffuse1;
	  _diffuse = diffuse(ray_origin, ray_direction, center1, t1, light_position);
          _diffuse = max(_diffuse, 0.); 
          _diffuse1 = diffuse(ray_origin, ray_direction, center1, t1, light_position1);
          _diffuse1 = max(_diffuse1, 0.); 
          _diffuse = _diffuse + _diffuse1;
          double _specular;
          _specular = specular(ray_origin, ray_direction, center1, t1, light_position, origin);
          _specular = pow(max(_specular, 0.),1000);
          R(i,j) =0*( _diffuse + _ambient + _specular);
          G(i,j) =0*( _diffuse + _ambient + _specular);
          B(i,j) =1*(_diffuse + _ambient);
          // Disable the alpha mask for this pixel
          A(i,j) = 1;


          }

          if(index == 2)
          {
          double _diffuse, _diffuse1;
	  _diffuse = diffuse(ray_origin, ray_direction, center2, t2, light_position);
          _diffuse = max(_diffuse, 0.); 
          _diffuse1 = diffuse(ray_origin, ray_direction, center2, t2, light_position1);
          _diffuse1 = max(_diffuse1, 0.); 
         _diffuse = _diffuse + _diffuse1;
          double _specular, _specular1;
          _specular = specular(ray_origin, ray_direction, center2, t2, light_position, origin);
          _specular = pow(max(_specular, 0.),10);
          _specular1 = specular(ray_origin, ray_direction, center2, t2, light_position1, origin);
          _specular1 = pow(max(_specular1, 0.),10);
          _specular = _specular + _specular1;
         // cout << "  " << _specular;
          R(i,j) =0*( _diffuse + _ambient + _specular);
          G(i,j) =1*( _diffuse + _ambient + _specular);
         // cout <<  G(i,j) << " ";
          B(i,j) =0*(_diffuse + _ambient + _specular);

           // Disable the alpha mask for this pixel
           A(i,j) = 1;}
         
          if(index == 3)
          {
           double _diffuse, _diffuse1;
	  _diffuse = diffuse(ray_origin, ray_direction, center3, t3, light_position);
          _diffuse = max(_diffuse, 0.); 
          _diffuse1 = diffuse(ray_origin, ray_direction, center3, t3, light_position1);
          _diffuse1 = max(_diffuse1, 0.);
          _diffuse = _diffuse + _diffuse1; 
          double _specular, _specular1;
          _specular = specular(ray_origin, ray_direction, center3, t3, light_position, origin);
          _specular = pow(max(_specular, 0.),50);
          _specular1 = specular(ray_origin, ray_direction, center3, t3, light_position1, origin);
          _specular1 = pow(max(_specular1, 0.),50);
          _specular = _specular + _specular1;

          R(i,j) =1*( _diffuse + _specular + _ambient);
          G(i,j) =0*( _diffuse + _specular + _ambient);
          B(i,j) =0*(_diffuse + _specular + _ambient);


          // Disable the alpha mask for this pixel
          A(i,j) = 1;}
                              
      }
           
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);

}

vector<string> Split(string str, string delim)
{
   vector<string> res;
   if("" == str) return res;
   string strs = str + delim;
   int pos;
   int size = strs.size();

   for(int i = 0; i < size; i++)
   {
	pos = strs.find(delim, i);
        if(pos < size)
        {
		string s = strs.substr(i, pos-i);
                res.push_back(s);
                i = pos + delim.size() - 1;
	}
   }
   return res;

}

void read_off1(string path, MatrixXd& V, MatrixXd& F)
{
	string str;
	char c;
	ifstream infile;
	infile.open(path);
	int ver_num;
	int line_num;
	int i = 0;


	while(i < 2)
	{
	getline(infile, str);
	if(i == 1)
	{
    	vector<string> vec = Split(str, " ");
    	ver_num = stoi(vec[0]);
        //ver_num = 3;
    	line_num = stoi(vec[1]) + stoi(vec[0]) + 2; 
    
	}

	i++;
	}

	V = MatrixXd::Zero(ver_num, 3);
	F = MatrixXd::Zero(line_num-ver_num-2, 3);
       
	while(i < ver_num + 2)
	{
	getline(infile, str);
	vector<string> vec = Split(str, " ");
	//for(int m = 0; m < 3; m++)
	//{
  	 //V(i-2, m) = stof(vec[m])*10;
         //if(m == 2)
         //V(i-2, m) =  V(i-2, m) + 0.4;
         //if(m == 1)
         //V(i-2, m) =  V(i-2, m) - 0.9;
        // cout <<  "m: " << m<< " n: "<< i-2<< " V: "<<V(m, i-2) << " ";
	//}
          V(i-2, 0) = -(stof(vec[2])*3.5 - 0.4);
          V(i-2, 1) = stof(vec[0])*3.5 - 0.5;
          V(i-2, 2) = stof(vec[1])*3.5 + 0.2;
        // cout << endl;
	i++;

	}


	while(i < line_num)
	{
	getline(infile, str);
	vector<string> vec = Split(str, " ");
	for(int m = 0; m < 3; m++)
	{
 	  F(i-2-ver_num, m) = stoi(vec[m+1]);
	  //cout <<  "m: " << m<< " n: "<< i-2-ver_num<< " F: "<<F(m, i-2-ver_num) << " ";
	}
        
        
        //cout << "aaaaaa";
	i++;

	}
        

	infile.close();

}

void read_off2(string path, MatrixXd& V, MatrixXd& F)
{
	string str;
	char c;
	ifstream infile;
	infile.open(path);
	int ver_num;
	int line_num;
	int i = 0;


	while(i < 2)
	{
	getline(infile, str);
	if(i == 1)
	{
    	vector<string> vec = Split(str, " ");
    	ver_num = stoi(vec[0]);
        //ver_num = 3;
    	line_num = stoi(vec[1]) + stoi(vec[0]) + 2; 
    
	}

	i++;
	}

	V = MatrixXd::Zero(ver_num, 3);
	F = MatrixXd::Zero(line_num-ver_num-2, 3);
       
	while(i < ver_num + 2)
	{
	getline(infile, str);
	vector<string> vec = Split(str, " ");
	//for(int m = 0; m < 3; m++)
	//{
  	 //V(i-2, m) = stof(vec[m])*10;
         //if(m == 2)
         //V(i-2, m) =  V(i-2, m) + 0.4;
         //if(m == 1)
         //V(i-2, m) =  V(i-2, m) - 0.9;
        // cout <<  "m: " << m<< " n: "<< i-2<< " V: "<<V(m, i-2) << " ";
	//}
          V(i-2, 0) = stof(vec[0])*0.08 - 0.2;
          V(i-2, 1) = stof(vec[1])*0.08 + 0.5;
          V(i-2, 2) = stof(vec[2])*0.08 + 0.9;
        // cout << endl;
	i++;

	}


	while(i < line_num)
	{
	getline(infile, str);
	vector<string> vec = Split(str, " ");
	for(int m = 0; m < 3; m++)
	{
 	  F(i-2-ver_num, m) = stoi(vec[m+1]);
	}
        
        
	i++;

	}
        

	infile.close();

}
Vector3d tri_normal(Vector3d a, Vector3d b, Vector3d c)
{
    return (b-a).cross(c-a).normalized();
}

Vector3d intersect(Vector3d a, Vector3d b, Vector3d c, Vector3d ray_origin, Vector3d ray_direction)
{   //cout << a << endl;
    
    MatrixXd A1 = MatrixXd::Zero(3,3);
    Vector3d a1 = a - b;
    Vector3d b1 = a - c;
    Vector3d b2 = a - ray_origin;
    //Vector3d nor = a1.cross(b1);
    //if(nor.dot(ray_origin) != 0.0)
    //{
	for(int m = 0; m < 3; m++)
	{
    	A1(m,0) = a1(m);
    	A1(m,1) = b1(m);
    	A1(m,2) = ray_direction(m);
	}
    
    return A1.partialPivLu().solve(b2);
   // return A1.fullPivHouseholderQr().solve(b2);
   // }
   // Vector3d parall(-1,-1,-1);
   // return parall;
  

}


double tri_diffuse(Vector3d ray_origin, Vector3d ray_direction, Vector3d ray_normal, double t, Vector3d light_position)
{
                   Vector3d ray_intersection = ray_origin + t * ray_direction;
                   
                   return (light_position-ray_intersection).normalized().dot(ray_normal);
}

double tri_specular(Vector3d ray_origin, Vector3d ray_direction, Vector3d ray_normal, double t, Vector3d light_position,  Vector3d origin)
{

		Vector3d ray_intersection = ray_origin + t * ray_direction;
		// Compute normal at the intersection point
		Vector3d v = (origin - ray_intersection).normalized();
                Vector3d l = (light_position - ray_intersection).normalized();
                Vector3d h = (v+l).normalized();
                return ray_normal.dot(h);
} 


double bound_r(MatrixXd V, Vector3d &bound_center)
{

    double max_x = -1.0;
    double max_y = -1.0;
    double max_z = -1.0;
    
    double min_x = 100.0;
    double min_y = 100.0;
    double min_z = 100.0;
  
    for(int s = 0 ; s < V.rows(); s++)  
    {
    if(V(s, 0) > max_x)
    max_x = V(s, 0);
    if(V(s, 0) < min_x)
    min_x = V(s, 0);

    if(V(s, 1) > max_y)
    max_y = V(s, 1);
    if(V(s, 1) < min_y)
    min_y = V(s, 1);

    if(V(s, 2) > max_z)
    max_z = V(s, 2);
    if(V(s, 2) < min_z)
    min_z = V(s, 2);

    bound_center(0) = (max_x + min_x)/2;
    bound_center(1) = (max_y + min_y)/2;
    bound_center(2) = (max_z + min_z)/2;
    }
    Vector3d max1(max_x, max_y, max_z);
    Vector3d min1(min_x, min_y, min_z);
    
    return ((max1-min1).squaredNorm())/2;

}


void part4()
{
    std::cout << "Part 4: Simple ray tracer, one bunny and one bumpy cube with orthographic projection" << std::endl;
    MatrixXd V1;
    MatrixXd F1;
    read_off1("../data/bunny.off", V1, F1);
   
    MatrixXd V2;
    MatrixXd F2;
    read_off2("../data/bumpy_cube.off", V2, F2);
    
    //bounding caculate
    Vector3d bound_center1;
    double r1 = bound_r(V1, bound_center1);
    Vector3d bound_center2;
    double r2 = bound_r(V2, bound_center2);


    const std::string filename("part4.png");
    MatrixXd R = MatrixXd::Zero(800,800); 
    MatrixXd G = MatrixXd::Zero(800,800);
    MatrixXd B = MatrixXd::Zero(800,800);// Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
      
    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(1.7,-1,2.15);
    Vector3d x_displacement(0,2.0/R.cols(),0);
    Vector3d y_displacement(0,0,-2.0/R.rows());

    
      
    // one light sources
    const Vector3d light_position(0, 0,1.4);
    const Vector3d light_position1(-1,-1,1.4);
   
    int t = 0;
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {   
           

	    // Prepare the ray
	    
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(-1,0,-0.15);

            //bounding test
            double tb1 = deta_ans(ray_origin, ray_direction, bound_center1, r1);
            double tb2 = deta_ans(ray_origin, ray_direction, bound_center2, r2);
           
            double tnear = 1000;
            int face = 0; 
            int index = 0;
            if(tb1 > 0.0 && tb1 < 100.0)
            
            {	
		
       	    	//compute the normal of all triangles
            	for(int k = 0; k < F1.rows(); k++)
           	{ 
            
           		Vector3d ans = intersect(V1.row(F1(k, 0)), V1.row(F1(k, 1)), V1.row(F1(k, 2)), ray_origin, ray_direction);
            

           		double u1 = ans(0);
            		double v1 = ans(1);
            		double t1 = ans(2);
            
           	 	if(t1 < tnear && t1 > 0 && u1 >= 0.0 && v1 >= 0.0 && u1 + v1 <= 1.0)
            		{
             		//cout<< " tnear: "<< tnear;
	    	 	tnear = t1;
             		face = k;
            		}

           
           	}

               
          }
          if(tb2 > 0.0 && tb2 < 100.0)
          {
            	 for(int k = 0; k < F2.rows(); k++)
          	 { 

	    	Vector3d ans = intersect(V2.row(F2(k, 0)), V2.row(F2(k, 1)), V2.row(F2(k, 2)), ray_origin, ray_direction);
           

            	double u2 = ans(0);
            	double v2 = ans(1);
            	double t2 = ans(2);
            	if(t2 < tnear && t2 > 0 && u2 > 0 && v2 > 0 && u2 + v2 <= 1)
            	{
	    	tnear = t2;
             	face = k;
             	index = 1;
            	}

           	}
          

          
          }

                // Final Shading model
           	if(tnear < 1000 && index == 0)
          	{
         	 Vector3d tnormal = tri_normal(V1.row(F1(face, 0)), V1.row(F1(face, 2)), V1.row(F1(face, 1)));
	  	double _diffuse, _diffuse1;
	  	_diffuse = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position);
          	_diffuse = max(_diffuse, 0.);
                _diffuse1 = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;
         
          	double _specular,_specular1 ;
          	_specular = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position, origin);
          	_specular = pow(max(_specular, 0.),1000);
		_specular1 = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position1, origin);
          	_specular1 = pow(max(_specular1, 0.),1000);
                _specular = _specular + _specular1;

	  	R(i,j) =0*( _diffuse + _ambient + _specular);
          	G(i,j) =0*( _diffuse + _ambient + _specular);
         	B(i,j) =1*(_diffuse + _ambient + _specular);
         	 // Disable the alpha mask for this pixel
          	A(i,j) = 1;

        	 }
                
                if(index == 1)
          	{
          	Vector3d tnormal = tri_normal(V2.row(F2(face, 0)), V2.row(F2(face, 1)), V2.row(F2(face, 2)));
	  	double _diffuse, _diffuse1;
	  	_diffuse = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position);
          	_diffuse = max(_diffuse, 0.);
		_diffuse1 = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;

          	double _specular, _specular1;
          	_specular = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position, origin);
          	_specular = pow(max(_specular, 0.),1000);
		_specular1 = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position1, origin);
          	_specular1 = pow(max(_specular1, 0.),1000);
                _specular = _specular + _specular1;

	  	R(i,j) =0.5*( _diffuse + _ambient + _specular);
          	G(i,j) =0*( _diffuse + _ambient + _specular);
          	B(i,j) =0.5*(_diffuse + _ambient + _specular);
          	// Disable the alpha mask for this pixel
          	A(i,j) = 1;
          	}
         
      	}
           
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);


}


bool intersect_object (Vector3d indot, Vector3d light_position, MatrixXd V, MatrixXd F)
{
		
		int intersecting = 0;
            	int k = 0;
                while(k < F.rows() && intersecting == 0)
           	{ 
			
           	Vector3d ans = intersect(V.row(F(k, 0)), V.row(F(k, 1)), V.row(F(k, 2)), indot, light_position - indot);
                k = k + 1;
        	double u1 = ans(0);
            	double v1 = ans(1);
            	double t1 = ans(2);
            	//intersecting with object and is shadow
           	if(t1 > 0 && u1 >= 0.0 && v1 >= 0.0 && u1 + v1 <= 1.0 && t1 < 1)
            	{
                   intersecting = 1;
     		   return true;
 		}
                }
                if(intersecting == 0)
                return false;
                
                return false;                    
}



void part5()
{
 std::cout << "Part 5: Simple ray tracer, one bunny and one cube with orthographic projection" << std::endl;
    MatrixXd V1;
    MatrixXd F1;
    read_off1("../data/bunny.off", V1, F1);
    
    MatrixXd V2;
    MatrixXd F2;
    read_off2("../data/bumpy_cube.off", V2, F2);
    
    //bounding caculate
    Vector3d bound_center1;
    double r1 = bound_r(V1, bound_center1);
    Vector3d bound_center2;
    double r2 = bound_r(V2, bound_center2);


    const std::string filename("part5.png");
    MatrixXd R = MatrixXd::Zero(400,400); 
    MatrixXd G = MatrixXd::Zero(400,400);
    MatrixXd B = MatrixXd::Zero(400,400);// Store the color
    MatrixXd A = MatrixXd::Zero(400,400); // Store the alpha mask
      
    // The camera is orthographic
    Vector3d origin(1.7,-1,2.15);
    Vector3d x_displacement(0,2.0/R.cols(),0);
    Vector3d y_displacement(0,0,-2.0/R.rows());

    //ground
    Vector3d planeA(6, 0, -0.6);
    Vector3d planeB(-5, 5, 1.2);
    Vector3d planeC(-5, -5,1.2);

      
    // one light sources
    const Vector3d light_position(0, 0,1.4);
    const Vector3d light_position1(-1,-1,1.4);

    int t = 0;
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {   
           

	    // Prepare the ray
	    
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(-1,0,-0.15);

            
            Vector3d ans = intersect(planeA, planeB, planeC, ray_origin, ray_direction);
            double t4 = ans(2);
            double u4 = ans(0);
            double v4 = ans(1);
            
            Vector3d indot = ray_origin + t4 * ray_direction;
            indot = indot + 1.0e-5 * (light_position - indot);
	    Vector3d indot1 = indot + 1.0e-5 * (light_position1 - indot);

            //bounding test
            double tb1 = deta_ans(ray_origin, ray_direction, bound_center1, r1);
            double tb2 = deta_ans(ray_origin, ray_direction, bound_center2, r2);     

            
                      
            double tnear = 1000;
            int face = 0; 
            int index = 0;
            if(tb1 > 0.0 && tb1 < 100.0)
            
            {	
		
       	    	//compute the normal of all triangles
            	for(int k = 0; k < F1.rows(); k++)
           	{ 
            
           		Vector3d ans = intersect(V1.row(F1(k, 0)), V1.row(F1(k, 1)), V1.row(F1(k, 2)), ray_origin, ray_direction);
            

           		double u1 = ans(0);
            		double v1 = ans(1);
            		double t1 = ans(2);
            
           	 	if(t1 < tnear && t1 > 0 && u1 >= 0.0 && v1 >= 0.0 && u1 + v1 <= 1.0)
            		{
             		//cout<< " tnear: "<< tnear;
	    	 	tnear = t1;
             		face = k;
            		}

           
           	}

               
          }
          if(tb2> 0.0 && tb2 < 100.0)
          {
            	 for(int k = 0; k < F2.rows(); k++)
          	 { 

	    	Vector3d ans = intersect(V2.row(F2(k, 0)), V2.row(F2(k, 1)), V2.row(F2(k, 2)), ray_origin, ray_direction);
           

            	double u2 = ans(0);
            	double v2 = ans(1);
            	double t2 = ans(2);
            	if(t2 < tnear && t2 > 0 && u2 > 0 && v2 > 0 && u2 + v2 <= 1)
            	{
	    	tnear = t2;
             	face = k;
             	index = 1;
            	}

           	}
          

          
          }
          if(tnear > 100)
          {
     
          bool shadow1 = intersect_object (indot, light_position, V1, F1);
          bool shadow2 = intersect_object (indot, light_position, V2, F2);
          bool shadow3 = intersect_object (indot1, light_position1, V1, F1);
          bool shadow4 = intersect_object (indot1, light_position1, V2, F2);
             	  if(shadow1 == false && shadow2 == false && shadow3 == false && shadow4 == false)	
                   {
                        R(i,j) =2*_ambient;
          		G(i,j) =2*_ambient;
          		B(i,j) =2*_ambient;
                   }
                       
                   else
                   {
                        R(i,j) =1*_ambient;
          	        G(i,j) =1*_ambient;
          	        B(i,j) =1*_ambient;
                        
                   }
           // Disable the alpha mask for this pixel
           A(i,j) = 1;
          }

                // Final Shading model
           	if(tnear < 1000 && index == 0)
          	{
         	 Vector3d tnormal = tri_normal(V1.row(F1(face, 0)), V1.row(F1(face, 2)), V1.row(F1(face, 1)));
	  	double _diffuse, _diffuse1;
	  	_diffuse = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position);
          	_diffuse = max(_diffuse, 0.);
                _diffuse1 = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;
         
          	double _specular,_specular1 ;
          	_specular = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position, origin);
          	_specular = pow(max(_specular, 0.),1000);
		_specular1 = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position1, origin);
          	_specular1 = pow(max(_specular1, 0.),1000);
                _specular = _specular + _specular1;

	  	R(i,j) =0*( _diffuse + _ambient + _specular);
          	G(i,j) =0*( _diffuse + _ambient + _specular);
         	B(i,j) =1*(_diffuse + _ambient + _specular);
         	 // Disable the alpha mask for this pixel
          	A(i,j) = 1;
        	 }
                
                if(index == 1)
          	{
          	Vector3d tnormal = tri_normal(V2.row(F2(face, 0)), V2.row(F2(face, 1)), V2.row(F2(face, 2)));
	  	double _diffuse, _diffuse1;
	  	_diffuse = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position);
          	_diffuse = max(_diffuse, 0.);
		_diffuse1 = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;

          	double _specular, _specular1;
          	_specular = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position, origin);
          	_specular = pow(max(_specular, 0.),1000);
		_specular1 = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position1, origin);
          	_specular1 = pow(max(_specular1, 0.),1000);
                _specular = _specular + _specular1;

	  	R(i,j) =0.5*( _diffuse + _ambient + _specular);
          	G(i,j) =0*( _diffuse + _ambient + _specular);
          	B(i,j) =0.5*(_diffuse + _ambient + _specular);
          	// Disable the alpha mask for this pixel
          	A(i,j) = 1;
          	}
         
      	}
           
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);





}






void part6()
{

std::cout << "Part 6: Simple ray tracer, one bunny and one cube with orthographic projection" << std::endl;
    MatrixXd V1;
    MatrixXd F1;
    read_off1("../data/bunny.off", V1, F1);
   
    MatrixXd V2;
    MatrixXd F2;
    read_off2("../data/bumpy_cube.off", V2, F2);
    
    //bounding caculate
    Vector3d bound_center1;
    double r1 = bound_r(V1, bound_center1);
    Vector3d bound_center2;
    double r2 = bound_r(V2, bound_center2);


    const std::string filename("part6.png");
    MatrixXd R = MatrixXd::Zero(400,400); 
    MatrixXd G = MatrixXd::Zero(400,400);
    MatrixXd B = MatrixXd::Zero(400,400);// Store the color
    MatrixXd A = MatrixXd::Zero(400,400); // Store the alpha mask
      
    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(1.7,-1,2.15);
    Vector3d x_displacement(0,2.0/R.cols(),0);
    Vector3d y_displacement(0,0,-2.0/R.rows());

    //ground
    Vector3d planeA(6, 0, -0.6);
    Vector3d planeB(-5, 5, 1.2);
    Vector3d planeC(-5, -5,1.2);

    //center of a circle
    Vector3d center2(-0.8,0.0,0.85);
    Vector3d center3(-1.4,-0.4,0.9);


    // one light sources
    const Vector3d light_position(0, 0,1.4);
    const Vector3d light_position1(-1,-1,1.4);
   
    int t = 0;
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {   
            const double sphere_radius2 = 0.2;
            const double sphere_radius3 = 0.4;

	    // Prepare the ray
	    
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(-1,0,-0.15);
          
            double sphere_t2 = deta_ans(ray_origin, ray_direction, center2, sphere_radius2);
            double sphere_t3 = deta_ans(ray_origin, ray_direction, center3, sphere_radius3);


            
            Vector3d ans = intersect(planeA, planeB, planeC, ray_origin, ray_direction);
            double t4 = ans(2);
            double u4 = ans(0);
            double v4 = ans(1);
            
            Vector3d indot = ray_origin + t4 * ray_direction;
            indot = indot + 1.0e-5 * (light_position - indot);
            Vector3d indot1 = indot + 1.0e-5 * (light_position1 - indot);
            Vector3d normal_dot = tri_normal(planeA, planeB, planeC);
            Vector3d ref = ray_direction - 2 * ray_direction.dot(normal_dot) * normal_dot;
           
            
            
            //bounding test
            double tb1 = deta_ans(ray_origin, ray_direction, bound_center1, r1);
            double tb2 = deta_ans(ray_origin, ray_direction, bound_center2, r2);     

            
            double tnear = 1000;
            int face = 0; 
            int index = 0;
            if(tb1 > 0.0 && tb1 < 100.0)
            
            {	
		
       	    	//compute the normal of all triangles
            	for(int k = 0; k < F1.rows(); k++)
           	{ 
            
           		Vector3d ans = intersect(V1.row(F1(k, 0)), V1.row(F1(k, 1)), V1.row(F1(k, 2)), ray_origin, ray_direction);
            

           		double u1 = ans(0);
            		double v1 = ans(1);
            		double t1 = ans(2);
            
           	 	if(t1 < tnear && t1 > 0 && u1 >= 0.0 && v1 >= 0.0 && u1 + v1 <= 1.0)
            		{
             		//cout<< " tnear: "<< tnear;
	    	 	tnear = t1;
             		face = k;
            		}

           
           	}

               
          }
          if(tb2> 0.0 && tb2 < 100.0)
          {
            	 for(int k = 0; k < F2.rows(); k++)
          	 { 

	    	Vector3d ans = intersect(V2.row(F2(k, 0)), V2.row(F2(k, 1)), V2.row(F2(k, 2)), ray_origin, ray_direction);
           

            	double u2 = ans(0);
            	double v2 = ans(1);
            	double t2 = ans(2);
            	if(t2 < tnear && t2 > 0 && u2 > 0 && v2 > 0 && u2 + v2 <= 1)
            	{
	    	tnear = t2;
             	face = k;
             	index = 1;
            	}

           	}
          }
          if(sphere_t2 < 100.0)
          {
            if(sphere_t2 < tnear)
		{
		tnear = sphere_t2;
                index = 2;
		}



	  }
          if(sphere_t3 < 100.0)
          {
            if(sphere_t3 < tnear)
		{
		tnear = sphere_t3;
                index = 3;
		}
          
          }

          if(tnear > 100 && sphere_t2 > 99.99 && sphere_t3 > 99.99)
          {

		double tnear1 = 100;
                int face1 = 0; 
                int index1 = 0;
                
                double tr = deta_ans(indot, ref, bound_center1, r1);
                double tr1 = deta_ans(indot, ref, bound_center2, r2);
                double tr2 = deta_ans(indot, ref, center2, sphere_radius2);
                double tr3 = deta_ans(indot, ref, center3, sphere_radius3);

                //if the refection ray has intersect with the bound, test all the face of the object
                if(tr > 0.0 && tr < 100.0)
            
            	{	
		
       	    	//compute the normal of all triangles
            	for(int k = 0; k < F1.rows(); k++)
           	{ 
            
           		Vector3d ans = intersect(V1.row(F1(k, 0)), V1.row(F1(k, 2)), V1.row(F1(k, 1)), indot, ref);
            

           		double u1 = ans(0);
            		double v1 = ans(1);
            		double t1 = ans(2);
            
           	 	if(t1 < tnear1 && t1 > 0 && u1 >= 0.0 && v1 >= 0.0 && u1 + v1 <= 1.0)
            		{
	    	 	tnear1 = t1;
             		face1 = k;
                        index1 = 1;
            		}

           
           	}

               
          	}
         	 if(tr1> 0.0 && tr1 < 100.0)
          	{
            	 for(int k = 0; k < F2.rows(); k++)
          	 { 

	    	Vector3d ans = intersect(V2.row(F2(k, 0)), V2.row(F2(k, 1)), V2.row(F2(k, 2)), indot, ref);
           

            	double u2 = ans(0);
            	double v2 = ans(1);
            	double t2 = ans(2);
            	if(t2 < tnear1 && t2 > 0 && u2 > 0 && v2 > 0 && u2 + v2 <= 1)
            	{
	    	tnear1 = t2;
             	face1 = k;
             	index1 = 2;
            	}
                
           	}
          	
          	}
                
		if(tr2 < tnear1)
		{
		tnear1 = tr2;
                index1 = 3;
		}
                if(tr3 < tnear1)
		{
		tnear1 = tr2;
                index1 = 4;
		}



                // shading the reflection
           	if(index1 == 1)
          	{
                  //cout << " face1: " <<face1;
                 //cout << " index1: " << index1;
         	 Vector3d tnormal = tri_normal(V1.row(F1(face1, 0)), V1.row(F1(face1, 2)), V1.row(F1(face1, 1)));
	  	double _diffuse, _diffuse1;
	  	_diffuse = tri_diffuse(indot, ref, tnormal, tnear1, light_position);
          	_diffuse = max(_diffuse, 0.); 
                _diffuse1 = tri_diffuse(indot1, ref, tnormal, tnear1, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;
       
          	double _specular,_specular1;
          	_specular = tri_specular(indot, ref, tnormal, tnear1, light_position, indot);
          	_specular = pow(max(_specular, 0.),200);
		_specular1 = tri_specular(indot1, ref, tnormal, tnear1, light_position1, indot);
          	_specular1 = pow(max(_specular1, 0.),200);
                _specular = _specular + _specular1;
	
         	B(i,j) =0.5*1*(_diffuse + _ambient + _specular);
         	
        	 }
                
                if(index1 == 2)
          	{
                 //cout << " face1: " <<face1;
                 //cout << " index1 : " << index1;
          	Vector3d tnormal = tri_normal(V2.row(F2(face1, 0)), V2.row(F2(face1, 1)), V2.row(F2(face1, 2)));
	  	double _diffuse,  _diffuse1;
	  	_diffuse = tri_diffuse(indot, ref, tnormal, tnear1, light_position);
          	_diffuse = max(_diffuse, 0.); 
                _diffuse1 = tri_diffuse(indot1, ref, tnormal, tnear1, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;

          	double _specular, _specular1;
          	_specular = tri_specular(indot, ref, tnormal, tnear1, light_position, indot);
          	_specular = pow(max(_specular, 0.),200);
   		_specular1 = tri_specular(indot1, ref, tnormal, tnear1, light_position1, indot);
          	_specular1 = pow(max(_specular1, 0.),200);
                _specular = _specular + _specular1;
 
	  	R(i,j) =0.5*1*( _diffuse + _ambient + _specular);
          	B(i,j) =0.5*1*(_diffuse + _ambient + _specular);
          	}


                if(index1 == 3)
		{
		double _diffuse, _diffuse1;
	  	_diffuse = diffuse(indot, ref, center2, tnear1, light_position);
          	_diffuse = max(_diffuse, 0.); 
          	_diffuse1 = diffuse(indot1, ref, center2, tnear1, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;
     
                G(i,j) =0.5*1*( _diffuse + _ambient);
		}


		if(index1 == 4)
		{
		double _diffuse, _diffuse1;
	  	_diffuse = diffuse(indot, ref, center3, tnear1, light_position);
          	_diffuse = max(_diffuse, 0.); 
          	_diffuse1 = diffuse(indot1, ref, center3, tnear1, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;
                double _specular, _specular1;
          	_specular = 1.1*specular(indot, ref, center3, tnear1, light_position, indot);
          	_specular = pow(max(_specular, 0.),200);
          	_specular1 = 1.1*specular(indot1, ref, center3, tnear1, light_position1, indot1);
         	_specular1 = pow(max(_specular1, 0.),200);
          	_specular = _specular + _specular1;
     
                R(i,j) =0.6 * 1 * ( _diffuse + _ambient + _specular);
		}
          
	  double t6 = deta_ans(indot, light_position - indot , center2, sphere_radius2);
          double t7 = deta_ans(indot, light_position - indot , center3, sphere_radius3);
          double t8 = deta_ans(indot1, light_position1 - indot1 , center2, sphere_radius2);
          double t9 = deta_ans(indot1, light_position1 - indot1 , center3, sphere_radius3);
     
          bool shadow1 = intersect_object (indot, light_position, V1, F1);
          bool shadow2 = intersect_object (indot, light_position, V2, F2);
	  bool shadow3 = intersect_object (indot1, light_position1, V1, F1);
          bool shadow4 = intersect_object (indot1, light_position1, V2, F2);
             	  if(shadow1 == false && shadow2 == false && shadow3 == false && shadow4 == false && t6 >= 1.0 && t7 >= 1.0 && t8 >= 1.0 && t9 >= 1.0)	
                   {
                        R(i,j) +=2*_ambient;
          		G(i,j) +=2*_ambient;
          		B(i,j) +=2*_ambient;
                   }
                       
                   else
                   {
                        R(i,j) +=1*_ambient;
          	        G(i,j) +=1*_ambient;
          	        B(i,j) +=1*_ambient;
                        
                   }
           // Disable the alpha mask for this pixel
           A(i,j) = 1;
          }

                // Final Shading model
           	if(tnear < 1000 && index == 0)
          	{
         	 Vector3d tnormal = tri_normal(V1.row(F1(face, 0)), V1.row(F1(face, 2)), V1.row(F1(face, 1)));
	  	double _diffuse, _diffuse1;
	  	_diffuse = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position);
          	_diffuse = max(_diffuse, 0.);
                _diffuse1 = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;
         
          	double _specular,_specular1 ;
          	_specular = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position, origin);
          	_specular = pow(max(_specular, 0.),200);
		_specular1 = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position1, origin);
          	_specular1 = pow(max(_specular1, 0.),200);
                _specular = _specular + _specular1;

	  	R(i,j) =0*( _diffuse + _ambient + _specular);
          	G(i,j) =0*( _diffuse + _ambient + _specular);
         	B(i,j) =1*(_diffuse + _ambient + _specular);
         	 // Disable the alpha mask for this pixel
          	A(i,j) = 1;
        	 }
                
                if(index == 1)
          	{
          	Vector3d tnormal = tri_normal(V2.row(F2(face, 0)), V2.row(F2(face, 1)), V2.row(F2(face, 2)));
	  	double _diffuse, _diffuse1;
	  	_diffuse = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position);
          	_diffuse = max(_diffuse, 0.);
		_diffuse1 = tri_diffuse(ray_origin, ray_direction, tnormal, tnear, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
                _diffuse = _diffuse + _diffuse1;

          	double _specular, _specular1;
          	_specular = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position, origin);
          	_specular = pow(max(_specular, 0.),200);
		_specular1 = tri_specular(ray_origin, ray_direction, tnormal, tnear, light_position1, origin);
          	_specular1 = pow(max(_specular1, 0.),200);
                _specular = _specular + _specular1;

	  	R(i,j) =0.5*( _diffuse + _ambient + _specular);
          	G(i,j) =0*( _diffuse + _ambient + _specular);
          	B(i,j) =0.5*(_diffuse + _ambient + _specular);
          	// Disable the alpha mask for this pixel
          	A(i,j) = 1;
          	}

                if(index ==2)
          	{
          	double _diffuse, _diffuse1;
	  	_diffuse = diffuse(ray_origin, ray_direction, center2, tnear, light_position);
          	_diffuse = max(_diffuse, 0.); 
          	_diffuse1 = diffuse(ray_origin, ray_direction, center2, tnear, light_position1);
          	_diffuse1 = max(_diffuse1, 0.); 
         	_diffuse = _diffuse + _diffuse1;
          	double _specular, _specular1;
          	_specular = specular(ray_origin, ray_direction, center2, tnear, light_position, ray_origin);
          	_specular = pow(max(_specular, 0.),100);
          	_specular1 = specular(ray_origin, ray_direction, center2, tnear, light_position1, ray_origin);
          	_specular1 = pow(max(_specular1, 0.),100);
          	_specular = _specular + _specular1;

         	 R(i,j) =0*( _diffuse + _ambient + _specular);
        	 G(i,j) =1*( _diffuse + _ambient +  _specular);
                 B(i,j) =0*(_diffuse + _ambient + _specular);

           	// Disable the alpha mask for this pixel
           	A(i,j) = 1;
                }


                if(index == 3)
		{
		 double _diffuse, _diffuse1;
	  	_diffuse = diffuse(ray_origin, ray_direction, center3, tnear, light_position);
          	_diffuse = max(_diffuse, 0.); 
          	_diffuse1 = diffuse(ray_origin, ray_direction, center3, tnear, light_position1);
          	_diffuse1 = max(_diffuse1, 0.);
          	_diffuse = _diffuse + _diffuse1; 
          	double _specular, _specular1;
          	_specular = 1.1*specular(ray_origin, ray_direction, center3, tnear, light_position, ray_origin);
          	_specular = pow(max(_specular, 0.),500);
          	_specular1 = 1.1*specular(ray_origin, ray_direction, center3, tnear, light_position1, ray_origin);
          	_specular1 = pow(max(_specular1, 0.),500);
          	_specular = _specular + _specular1;

         	R(i,j) =1*( _diffuse + _ambient+ _specular);
          	G(i,j) =0*( _diffuse + _specular + _ambient);
          	B(i,j) =0*(_diffuse + _specular + _ambient);


          	// Disable the alpha mask for this pixel
          	A(i,j) = 1;

		}
         
      	}
           
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);





}



int main()
{
 
    part1();
    part2();
    part3();
    part4();
    part5();
    part6();

    return 0;
}
