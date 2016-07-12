/* 
   bond_angle.c
   Usage: [binary] <filename.lammpstrj>
   Compilation requires C++11/C++0x (g++ bond_angle.c -std=c++0x)

   Calculates the orientation of bond angles for LAMMPS simulations 
   involving triatomic molecules (e.g. water). 
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>		/* find */

#define PRELIM_LINES 9 		/* Metadata lines in lammpstrj */

using namespace std;

typedef vector<double> vd;

void print_vector(vd vec);
void adjacency(vector<int> atom_id, vector<int> type_id,
	       vector<int> &adj0, vector<int> &adj1, vector<int> &adj2);
void get_metadata(string filename, int &numAtoms, int &timestep, int &numFrames);
void build_adj_list(vector<string> raw_file, int frame_count, int num_atoms,
		    vector<int> &atom_id, vector<int> &type_id,
		    vector<int> &adj0, vector<int> &adj1, vector<int> &adj2);

double norm(double x, double y, double z);
void calculate_angle(vector<string> raw_file, int frame_count, int num_atoms,
		     vector<int> adj0, vector<int> adj1, vector<int> adj2);

int main(int argc, char *argv[]) {

  /* No file specified */
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <filename>\n";
    return 1;
  }
  else {
    /* Check file extension */
    string filename = argv[1];
    string to_string(filename);
    int ext_pos = filename.find_last_of(".");
    string file_ext = filename.substr(ext_pos+1);
      
    if (file_ext != "lammpstrj") {
      cout << "Error: not a LAMMPS trajectory. \n";
      return 1;
    } else {
      /* The file extension is correct, open it */
      string line;
      ifstream myfile;
      int count = 0;

      /* Metadata from the Lammps trajectory */
      int num_atoms = 0;	/* number of atoms */
      int time_step = 0;      	/* interval between traj dump */
      int num_frames = 0;	/* number of frames in the traj */
      get_metadata(argv[1], num_atoms, time_step, num_frames);
      
      myfile.open(argv[1]);
      
      if(!myfile) {
	cout << "Error: unable to open trajectory. \n";
	return 1;
      } else {

	/* Save the raw data into memory */
	vector<string> raw_file; 
	if (myfile.is_open()) {
	  while(getline(myfile, line)) { raw_file.push_back(line); }
	  myfile.close();
	} else { cout << "Error: unable to open file. \n";}

	/* Populates an adjacency matrix from the first timestep */
	vector<int> atom_id_0, type_id_0;	   /* atom id and type id */
	vector<double> x, y, z;	/* coordinate of each atom */
	vector<int> adj0, adj1, adj2; 
	build_adj_list(raw_file, 0, num_atoms, atom_id_0, type_id_0, 
		       adj0, adj1, adj2);

	/* We'll use the adjacency matrix to track the molecules */
	for (int i=0; i<2500; i++) {
	  calculate_angle(raw_file, i+2500, num_atoms, adj0, adj1, adj2);
	  cout << i << endl;
	}
 
      }
    }    
  }
}

/* Read key data from the trajectory */
void get_metadata(string filename, int &numAtoms,
		  int &time_step, int &numFrames) {
  ifstream meta;
  string line;
  int count = 0;
  int ts_line;
  meta.open(filename);
  if (meta.is_open()) {
    while(getline(meta, line)) {
      if (count == 3) {
	numAtoms = stoi(line);
	ts_line = numAtoms + PRELIM_LINES + 1;
      }
      if (count == ts_line) {
	time_step = stoi(line);
      }
      count += 1; /* iterate line counter */
    }
  }
  meta.close();
  numFrames = count/(numAtoms+PRELIM_LINES);
  return;
}

void print_vector(vd vec) {
  for (int i=0; i<vec.size(); i++) { cout << vec[i] << endl;}
}


/* For a frame_count (out of num_frames), calculate the bond angles */
void build_adj_list(vector<string> raw_file, int frame_count, int num_atoms,
		    vector<int> &atom_id, vector<int> &type_id,
		    vector<int> &adj0, vector<int> &adj1, vector<int> &adj2) {

  int first_line = PRELIM_LINES + frame_count*(PRELIM_LINES + num_atoms);
  int last_line = (frame_count+1)*(PRELIM_LINES + num_atoms);
  string s;
  int count;

  /* Populate vectors with the coordinates of each atom */
  for (int i=first_line; i<last_line; i++) {
    count = 0;
    istringstream ss (raw_file[i]);
    while (getline(ss, s, ' ')) {
      if (count == 0) {atom_id.push_back(stoi(s));}
      if (count == 1) {type_id.push_back(stoi(s));}      
      count += 1;
    }
  }

  adjacency(atom_id, type_id, adj0, adj1, adj2);
  
  return; 
}

/* Calculate adjacency vector */
void adjacency(vector<int> atom_id, vector<int> type_id,
	       vector<int> &adj0, vector<int> &adj1, vector<int> &adj2) {

  int cline = 0;
  for (int i=0; i< atom_id.size(); i++) {
    cline = i*3;    
    if (type_id[cline] == 2) {	
      adj0.push_back(atom_id[cline]);    
      adj1.push_back(atom_id[cline+1]);
      adj2.push_back(atom_id[cline+2]);
    }
   }
  
  return;
}

/* Calculates the norm of a vector */
double norm(double x, double y, double z) {
  return sqrt(x*x + y*y + z*z);
}


void calculate_angle(vector<string> raw_file, int frame_count, int num_atoms,
		     vector<int> adj0, vector<int> adj1, vector<int> adj2) {
  
  int first_line = PRELIM_LINES + frame_count*(PRELIM_LINES + num_atoms);
  int last_line = (frame_count+1)*(PRELIM_LINES + num_atoms);
  string s;
  int count;

  vector<int> atom_id, type_id;
  vector<double> x, y, z;

  /* Populate vectors with the coordinates of each atom */
  for (int i=first_line; i<last_line; i++) {
    count = 0;
    istringstream ss (raw_file[i]);
    while (getline(ss, s, ' ')) {
      if (count == 0) {atom_id.push_back(stoi(s));}
      if (count == 1) {type_id.push_back(stoi(s));}      
      if (count == 2) {x.push_back(stod(s)); }
      if (count == 3) {y.push_back(stod(s)); }
      if (count == 4) {z.push_back(stod(s)); }      
      count += 1;
    }
  }

  /* Sort the vectors according to the order given at time step 0 */
  /* so we can track them */
  vector<double> xs, ys, zs;
  int o1, h1, h2;
  for (int i=0; i<adj0.size(); i++) {
    o1 = find(atom_id.begin(), atom_id.end(), adj0[i]) - atom_id.begin();
    h1 = find(atom_id.begin(), atom_id.end(), adj1[i]) - atom_id.begin();
    h2 = find(atom_id.begin(), atom_id.end(), adj2[i]) - atom_id.begin();    

    xs.push_back(x[o1]); xs.push_back(x[h1]); xs.push_back(x[h2]);    
    ys.push_back(y[o1]); ys.push_back(y[h1]); ys.push_back(y[h2]);
    zs.push_back(z[o1]); zs.push_back(z[h1]); zs.push_back(z[h2]);    
    
  }

  /* Generate basis vectors for the plane formed by each trio of atoms */
  int cur_line;
  vector<double> v1_x, v1_y, v1_z; 
  vector<double> v2_x, v2_y, v2_z; 	
  
  for (int i=0; i<xs.size()/3; i++) {
    cur_line = 3*i;
    v1_x.push_back(xs[cur_line] - xs[cur_line+1]);
    v1_y.push_back(ys[cur_line] - ys[cur_line+1]);
    v1_z.push_back(zs[cur_line] - zs[cur_line+1]);
    v2_x.push_back(xs[cur_line] - xs[cur_line+2]);
    v2_y.push_back(ys[cur_line] - ys[cur_line+2]);
    v2_z.push_back(zs[cur_line] - zs[cur_line+2]);    
  }

  /* Perform a cross product */
  vector<double> n_x, n_y, n_z;
  for (int i=0; i < v1_x.size(); i++) {
    n_x.push_back(v1_y[i]*v2_z[i] - v1_z[i]*v2_y[i]);
    n_y.push_back(v1_z[i]*v2_x[i] - v1_x[i]*v2_z[i]);
    n_z.push_back(v1_x[i]*v2_y[i] - v1_y[i]*v2_x[i]);        
  }
  
  vector<double> ref_vec = {0.0, 1.0, 0.0};
  
  /* Calculate the angle */
  vector<double> angle;
  double ang_arg, vec_norm; 
  for (int i=0; i<n_x.size(); i++) {
    ang_arg = n_x[i]*ref_vec[0] + n_y[i]*ref_vec[1] + n_z[i]*ref_vec[2];
    vec_norm = norm(n_x[i], n_y[i], n_z[i]);
    angle.push_back(acos(ang_arg/(vec_norm)));
  }

  /* Output angle data */
  ofstream outfile;
  string outfilename;
  outfilename.append("file_");
  outfilename.append(to_string(frame_count));
  outfilename.append(".txt");
  outfile.open(outfilename);
  for (int i=0; i<angle.size(); i++)
    outfile << angle[i] << endl;
  outfile.close();
  
  return; 
}
