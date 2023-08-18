#include "../include/util.h"
#include "../include/hcal.h"

namespace util {

  ////Console
  // returns string with current month day and year (date)
  string getDate(){
    time_t now = time(0);
    tm ltm = *localtime(&now);
  
    string yyyy = to_string(1900 + ltm.tm_year);
    string mm = to_string(1 + ltm.tm_mon);
    string dd = to_string(ltm.tm_mday);
    string date = mm + '_' + dd + '_' + yyyy;
  
    return date;
  }

  // reads a set of parameters into a vector from a txt file at <const_path>. Assumes endl delimiter
  void readParam( std::string const_path, vector<Double_t> &param ) {

    ifstream const_file; const_file.open( const_path );

    //Get all params
    if(const_file.is_open()){
  
      Int_t n1=0;
      Double_t d1;
      string readline;
    
      while( getline( const_file, readline ) ){
	if( readline.at(0) == '#' ) {
	  continue;
	}
	istringstream iss( readline );
	iss >> d1;     
	param.push_back(d1);
	n1++;
      }
      
    }else{
      cerr << "Error on [util::readParam] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
    const_file.close();
  }

  // Counts the number of calibration sets for a given calibration parameter
  Int_t countSets( std::string const_path,
		   std::string type ){
    
    Int_t Nsets = 0;

    ifstream const_file; const_file.open( const_path.c_str() );

    //Get all params at datatype
    if(const_file.is_open()){

      string readline;

      while( getline( const_file, readline ) ){
	
	TString Tline = (TString)readline;
	
	if( Tline.BeginsWith(type) )
	  Nsets++;

      }

      return Nsets;

    }else{
      cerr << "Error on [util::countSets] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
  } 
  

  // Reads database-style file for 288 (hcal::maxHCalChan) elements after a given moniker (type)
  // A return character must exist after the list of parameters or error will be thrown
  void readDB( std::string const_path, 
	       std::string timestamp, 
	       std::string type,
	       Double_t param[] ) {

    ifstream const_file; const_file.open( const_path.c_str() );
    
    //Get all params at datatype
    if(const_file.is_open()){

      Int_t n0=0;
      Double_t d0;
      std::string readline;
      bool read_param = false;
      bool found_tstamp = false;

      if( timestamp.compare("none")==0 )
	found_tstamp = true;

      while( getline( const_file, readline ) ){
      
	if( n0==hcal::maxHCalChan ){
	  const_file.close();
	  return;
	}
      
	TString Tline = (TString)readline;
  
	if( Tline.BeginsWith(timestamp) && !found_tstamp ){
	  found_tstamp = true;
	  continue;
	}

	if( Tline.BeginsWith(type) && !read_param && found_tstamp ){
	  read_param = true;
	  continue;
	}
      
	if( found_tstamp && read_param ){

	  istringstream iss( readline );
	  while( iss >> d0 ){
	    
	    param[n0] = d0;

	    n0++;
	  }
	}
      }
      
      cerr << "Error on [util::readDB] Did not find required number of elements in param[] at " << const_path << endl;
      throw;

    }else{
      cerr << "Error on [util::readDB] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
  }

  // Overload for tdc calibration constant (1 element, same line)
  void readDB( std::string const_path, 
	       std::string timestamp, 
	       std::string type,
	       Double_t &param ) {

    ifstream const_file; const_file.open( const_path.c_str() );
    
    //Get all params at datatype
    if(const_file.is_open()){

      Double_t d0;
      string readline;
      bool found_tstamp = false;

      if( timestamp.compare("none")==0 )
	found_tstamp = true;

      while( getline( const_file, readline ) ){

	TString Tline = (TString)readline;
	
	if( Tline.BeginsWith(timestamp) && !found_tstamp ){	  
	  found_tstamp = true;
	  continue;
	}

	if( found_tstamp && Tline.BeginsWith(type) ){
	  string temp;
	  bool val = false;
	  for( Int_t l=0; l<readline.length(); l++){
	    if( val==true )
	      temp+=readline[l];
	    if( readline[l] == '=' ) 
	      val=true;
	  }
	  param = stof( temp );
	  const_file.close();
	  return;
	}
      }
      cerr << "Error on [util::readDB] timestamp/type does not exist in file at " << const_path << endl;
      throw;

    }else{
      cerr << "Error on [util::readDB] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
  }


  // Overload for other constants (N elements, same line)
  void readDB( std::string const_path, 
	       std::string timestamp, 
	       std::string type,
	       vector<Double_t> &param ) {

    ifstream const_file; const_file.open( const_path.c_str() );
    
    //Get all params at datatype
    if(const_file.is_open()){

      Double_t d0;
      string readline;
      bool found_tstamp = false;

      if( timestamp.compare("none")==0 )
	found_tstamp = true;

      while( getline( const_file, readline ) ){

	TString Tline = (TString)readline;
	
	if( Tline.BeginsWith(timestamp) && !found_tstamp ){	  
	  found_tstamp = true;
	  continue;
	}

	if( found_tstamp && Tline.BeginsWith(type) ){
	  string temp;
	  bool val = false;
	  for( Int_t l=0; l<readline.length(); l++){
	    if( val==true )
	      temp+=readline[l];
	    if( readline[l] == '=' ) 
	      val=true;
	  }
	  istringstream iss( temp );
	  while( iss >> d0 )
	    param.push_back(d0);

	  const_file.close();
	  return;
	}
      }
      cerr << "Error on [util::readDB] timestamp/type does not exist in file at " << const_path << endl;
      throw;

    }else{
      cerr << "Error on [util::readDB] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
  }

  // compares two timestamps and returns the one that is later in time. Format must be similar
  void tsCompare( std::string tsA, 
  		  std::string tsB, 
  		  std::string &tsLate ){
    
    std::string tempA;
    std::string tempB;

    if( tsA.compare(tsB)==0 ){ //if they are the same, return tsA arb.
      tsLate=tsA;
      return;
    }else if(tsA.compare("none")==0){ //if they are not, check if A is none, then return the other
      tsLate=tsB;
      return;
    }else if(tsB.compare("none")==0){ //vice versa
      tsLate=tsA;
      return;
    }else{ //else, get new strings deleting "-------[ " preamble as it can be of different length
      bool Atoi = false;
      for( Int_t l=0; l<tsA.size(); l++ ){
	char a = tsA[l];
	if( Atoi )
	  tempA+=a;
	if( a==' ' )
	  Atoi = true;
      }
      bool Btoi = false;
      for( Int_t l=0; l<tsB.size(); l++ ){
	char b = tsB[l];
	if( Btoi )
	  tempB+=b;
	if( b==' ' )
	  Btoi = true;
      }
    
      Int_t l=0;
      while( tempA[l++] ){ //compare the two parsed strings and return newest
	
    	char a = tempA[l];
    	char b = tempB[l]; 

	if( (Int_t)a>(Int_t)b ){
	  tsLate = tsA;
	  return;
	}else if( (Int_t)a<(Int_t)b ){
	  tsLate = tsB;
	  return;
	}else
	  continue;
      }

    }

    cerr << "Erroron [util::tsCompare] String comparison failed" << endl;
    throw; //comparison should be exhaustive

  }
  
  //Reads run list .csv file
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   std::string experiment,     // experiment {gmn,gen,genrp,gep,etc.}
		   Int_t &nruns,               // No. of runs to analyze
		   Int_t sbsconf,              // SBS configuration
		   Int_t replay_pass,          // replay pass
		   Int_t verbose,              // verbosity
		   vector<calrun> &crun)       // Output: Vector of crun structs
  {
    // Define the name of the relevant run spreadsheet
    if ( experiment.compare("gmn")==0 && replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = Form("%s%sruns_pass%d.csv",runsheet_dir.c_str(),experiment.c_str(),replay_pass);

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    ifstream run_data; run_data.open(run_spreadsheet);
    string readline;
    if(run_data.is_open()){
      std::cout << "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;     
      string skip_header; getline(run_data, skip_header); // skipping column header
      while(getline(run_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;
	  temp.push_back(temptoken);
	}
	// add relevant info to calrun objects
	if (stoi(temp[0]) == sbsconf ) {
	  if (crun.size() >= nruns) break;
	  calrun temp_cr;
	  temp_cr.SetDataRunSheet(temp);
	  crun.push_back(temp_cr);
	}

	temp.clear();
      }
      //Update nruns with total no. of runs to analyze
      nruns = crun.size();
      if (verbose == 1) {
	std::cout << "First run info:" << std::endl << crun[0];
	std::cout << "Last run info:" << std::endl << crun[nruns-1];
      }
    }else{
      cerr << "Error on [util::ReadRunList], run spreadsheet doesn't exist" << endl;
      throw;
    }
    run_data.close();
  }

  //Read cut list .csv file. sbsconf + target + field MUST uniquely specify cuts loaded here
  void ReadCutList(std::string cutsheet_dir,  // Dir. path containing CSV files with cut info
		   std::string experiment,     // experiment {gmn,gen,genrp,gep,etc.}
		   Int_t sbsconfig,            // SBS configuration
		   Int_t calib_set,            // Calibration set
		   Int_t replay_pass,          // replay pass
		   std::string target,         // target
		   Int_t field,                // sbs field in percent
		   Int_t verbose,              // verbosity
		   vector<calcut> &ccut)       // Output: Vector of ccut structs
  {
    // Define the name of the relevant cut spreadsheet
    if ( experiment.compare("gmn")==0 && replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string cut_spreadsheet = Form("%s%scuts_pass%d.csv",cutsheet_dir.c_str(),experiment.c_str(),replay_pass);

    // Reading the spreadsheet
    ifstream cut_data; cut_data.open(cut_spreadsheet);
    string readline;
    bool max_pushback = false;
    if(cut_data.is_open()){     
      string skip_header; getline(cut_data, skip_header); // skipping column header
      while(getline(cut_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;
	  temp.push_back(temptoken);
	}
	// select cut based on config, target, and field
	if (stoi(temp[0]) == sbsconfig &&
	    stoi(temp[1]) == calib_set &&
	    temp[2].compare(target)==0 &&
	    stoi(temp[3]) == field ) {
	  if( max_pushback ) //throw error if more than one element added to cut array
	    throw "Error on [util::ReadCutList], config/target/field not uniquely specified in cut spreadsheet";
	  calcut temp_cc;
	  temp_cc.SetDataCutSheet(temp);
	  ccut.push_back(temp_cc);
	  max_pushback = true;
	}

	temp.clear();
      }
      //Update nruns with total no. of runs to analyze
      if (verbose == 1) {
	std::cout << "First cut info:" << std::endl << ccut[0];
	std::cout << "Last cut info:" << std::endl << ccut.back();
      }
    }else{
      cerr << "Error on [util::ReadCutList], cut spreadsheet doesn't exist" << endl;
      throw;
    }
    cut_data.close();
  }

  //Geometry/physics
  //sets hcal axes by kinematic (from hall coordinates)
  void sethcalaxes( Double_t sbstheta_rad,                           // SBS angle in radian 
		    vector<TVector3> &hcal_axes ) {
    TVector3 hcal_zaxis( sin(-sbstheta_rad), 0, cos(-sbstheta_rad) ); // Clock-wise rotation about Y axis
    TVector3 hcal_xaxis( 0, -1, 0 );                                  // -Y axis of Hall CoS = X axis of hcal CoS
    TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
    hcal_axes.push_back(hcal_xaxis);
    hcal_axes.push_back(hcal_yaxis);
    hcal_axes.push_back(hcal_zaxis);
  }

  //gets expected location of scattered nucleon assuming straight line projections from BB track
  void getxyhcalexpect( TVector3 vertex, TVector3 pNhat, TVector3 hcal_origin, 
			vector<TVector3> hcal_axes, vector<Double_t> &xyhcalexpect ) {
    // intersection of a ray with a plane in hcal coordinates
    Double_t sintersect = ( hcal_origin - vertex).Dot(hcal_axes[2] ) / ( pNhat.Dot(hcal_axes[2]) );
    // ray from Hall origin onto the face of hcal where the nucleon hit
    TVector3 hcal_intersect = vertex + sintersect*pNhat; 

    Double_t xexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[0] );
    Double_t yexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[1] );

    xyhcalexpect.push_back( xexpect_hcal );
    xyhcalexpect.push_back( yexpect_hcal );
  }
  
  //checks particle id with information from observed dx distributions.
  void checkPID( std::string target, Double_t dx0p, Double_t dx0n, Double_t dxsigp, 
		 Double_t dxsign, Double_t dx, Int_t &pid ){
    bool isproton = abs(dx-dx0p)<3*dxsigp;
    bool isneutron = abs(dx-dx0n)<3*dxsign;
    bool isambiguous = isproton && isneutron;
    
    if( isproton ) 
      pid = 1;
    else if( isneutron )
      pid = 2;
    else if( isambiguous ){
      if( target.compare("lh2")==0 || target.compare("LH2") ) //probably not useful since isproton comes first
	pid = 1;
      else
	pid = 0;
    }
    else
      pid = -1; //neither
  }
  
  //Fit functions
  //gaussian fit
  Double_t g_gfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
  }
 
  //gaussian fit with uniform offset
  Double_t g_gfit_bg(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    Double_t bg = par[3];
    return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.)) + bg;
  }
 
  //skewed gaussian fit
  Double_t g_sgfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    Double_t alpha = par[3];
    return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) );
  }

  //skewed gaussian fit with uniform offset
  Double_t g_sgfit_bg(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    Double_t alpha = par[3];
    Double_t bg = par[4];
    return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) ) + bg;
  }

  //expo fit
  Double_t g_expofit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t str = par[1];
    return exp(amp+str*x[0]);
  }

  //traditional timewalk fit
  Double_t g_tradtwfit(Double_t *x, Double_t *par){
    Double_t asymp = par[0];
    Double_t scale = par[1];
    Double_t expo = par[2];
    return asymp+scale/pow(x[0],expo);
  }

  //tuned exponential timewalk fit
  Double_t g_twfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t str = par[1];
    Double_t offset = par[2];
    return amp*exp(-str*x[0])+offset;
  } 

  //expo fit with offset
  Double_t g_scexpofit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t amp = par[1];
    Double_t str = par[2];
    return yint+exp(amp+str*x[0]);
  }

  //linear fit
  Double_t g_lfit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];;
    return yint+p1*x[0];
  }

  //3rd order poly fit
  Double_t g_p3fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3);
  }

  //4th order poly fit
  Double_t g_p4fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4);
  }

  //6th order poly fit
  Double_t g_p6fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6);
  }

} //::util
