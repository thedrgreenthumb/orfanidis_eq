#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "orfanidis_eq.h"

using namespace std;
using namespace orfanidis_eq;

//Include testing data from file:
//sample rate in Hz, number of bands, number of samples in testing vectors
const unsigned int cfg_data[3] =
{
#include "orfanidis_eq.tstdat"
};

#define READED_SAMPLE_RATE (cfg_data[0])
#define READED_NUMBER_OF_BANDS (cfg_data[1])
#define TEST_VECTORS_NUMBER_OF_SAMPLES (cfg_data[2])

//Include testing data from file:
//30 gain values for equalizers configuration
const double cfg_gains_data[30] = //TODO:Find some other way to read gains
{
#include "orfanidis_eq_gain.tstdat"
};

static const int class_conversions_min_max_db = eq_min_max_gain_db;

static const double default_max_gain = 1.0;
static const double default_min_gain = 0.0;
static const double default_unit_imp_amp = 1.0;

// ------------ test_orfanidis_eq ------------
//Provide auxiliary functions for EQ's testing
template <typename eq_type>
class test_orfanidis_eq
{
public:
    test_orfanidis_eq(){}
    ~test_orfanidis_eq(){}

	bool set_freq_grid(freq_grid& fg, unsigned int number_of_bands) {
        switch(number_of_bands)
        {
            case 5:
                fg.set_5_bands(bands_grid_center_freq_hz);
                break;
            case 10:
                fg.set_10_bands(bands_grid_center_freq_hz);
                break;
            case 20:
                fg.set_20_bands(bands_grid_center_freq_hz);
                break;
            case 30:
                fg.set_30_bands(bands_grid_center_freq_hz);
                break;
            default:
                cout << "Can not configue freq grid, provide 5,10,20 or 30." << 					"Provided : " << number_of_bands << endl;
                return true;
        }
        return false;
    }

	void set_unit_impulse(vector<eq_single_t>& vect) {
		vect[0] = default_unit_imp_amp;
	}
	
	void process_eq(eq_type& equalizer, vector<eq_single_t> &in, 
		vector<eq_single_t> &out) {
        for(unsigned int i = 0; i < in.size(); i++)
            equalizer.sbs_process(&in[i], &out[i]);
    }
	
	void save_cs_file(string file_name, vector<eq_single_t>& test_vector) {
        //Write to file
        ofstream data_file;
        data_file.open(file_name.c_str());
        if (!data_file.is_open())
            return;
        for(unsigned int j = 0; j < test_vector.size(); j++)
            data_file << test_vector[j] << ", ";
        data_file.close();     
    }
};


void print_input_configuration_data() {
	//Print input cfg data
    cout << "Data from cfg file : " << endl;
	cout << READED_SAMPLE_RATE << " " << 
		READED_NUMBER_OF_BANDS << " " << TEST_VECTORS_NUMBER_OF_SAMPLES << endl;
	
	cout << "Gains cfg data : " << endl;
	for(unsigned int i = 0; i < READED_NUMBER_OF_BANDS; i++)
		cout << cfg_gains_data[i] << ", ";
		
	cout << endl << "dB" << endl;
}

void class_conversions_test() {
	//Conversions class test
	conversions convs(class_conversions_min_max_db);
	cout << "1. conversions class: fast conversions test : " << endl;
	cout << "db -> lin -> fast_lin" << endl;
	eq_double_t test_min_max_db = class_conversions_min_max_db + 3;
	eq_single_t ref_db = 0;
	for(ref_db = -test_min_max_db; ref_db < test_min_max_db; ref_db+=6.4)
		 cout << ref_db << " " << convs.db_2_lin(ref_db) << " " <<  	
		 	convs.fast_db_2_lin(ref_db) << endl;

	cout << "lin -> db -> fast_db" << endl; 
	eq_single_t ref_lin = 0;	
 	for(ref_lin = 0.00316 /* -50 db */; ref_lin < 316 /* 50 db */; 
 		ref_lin*=convs.db_2_lin(6))
 		cout << ref_lin << " " << convs.lin_2_db(ref_lin) << " " << 
		 	convs.fast_lin_2_db(ref_lin) << endl;
	
	cout << "----" << endl;
}

void class_freq_grid_test() {
	cout << "2. freq_grid class: print freqs for 1/3 octave eq:" << endl;
	freq_grid fg;
	fg.set_30_bands();
	
	cout << "band freq : rounded band freq" << endl;
	for (unsigned int i = 0; i < fg.get_number_of_bands(); i++)
	cout << fg.get_freq(i) << " : " << fg.get_rounded_freq(i) << endl;
}

template <typename eq_type>
void class_eq_test(string extention) {
	test_orfanidis_eq<eq_type> test_eq;
	
	//Input configuration
	freq_grid fg;
    test_eq.set_freq_grid(fg, READED_NUMBER_OF_BANDS);
    
    //Input data vector
    vector<eq_single_t> in_vector(TEST_VECTORS_NUMBER_OF_SAMPLES, 0);
    test_eq.set_unit_impulse(in_vector);
    
    //Input gains vector
    vector<eq_single_t> gains(cfg_gains_data, 
    	cfg_gains_data + READED_NUMBER_OF_BANDS);    

    //Butterworth
    vector<eq_single_t> butterworth_out(TEST_VECTORS_NUMBER_OF_SAMPLES, 0);
	eq_type equalizer(fg, butterworth);
	equalizer.set_sample_rate(READED_SAMPLE_RATE);
	equalizer.change_gains_db(gains);
	test_eq.process_eq(equalizer, in_vector, butterworth_out);
	string butterworth_f_name("butterworth_");
	butterworth_f_name += extention;
	butterworth_f_name += ".tstdat";
    test_eq.save_cs_file(butterworth_f_name , butterworth_out);
    
    //Chebyshev1
    vector<eq_single_t> chebyshev1_out(TEST_VECTORS_NUMBER_OF_SAMPLES, 0);
    equalizer.set_eq(fg, chebyshev1);
	equalizer.set_sample_rate(READED_SAMPLE_RATE);
	equalizer.change_gains_db(gains);
	test_eq.process_eq(equalizer, in_vector, chebyshev1_out);
	string chebyshev1_f_name("chebyshev1_");
	chebyshev1_f_name += extention;
	chebyshev1_f_name += ".tstdat";
    test_eq.save_cs_file(chebyshev1_f_name, chebyshev1_out);
    
    //Chebyshev2
    vector<eq_single_t> chebyshev2_out(TEST_VECTORS_NUMBER_OF_SAMPLES, 0);
    equalizer.set_eq(fg, chebyshev2);
	equalizer.set_sample_rate(READED_SAMPLE_RATE);
	equalizer.change_gains_db(gains);
	test_eq.process_eq(equalizer, in_vector, chebyshev2_out);
	string chebyshev2_f_name("chebyshev2_");
	chebyshev2_f_name += extention;
	chebyshev2_f_name += ".tstdat";
    test_eq.save_cs_file(chebyshev2_f_name, chebyshev2_out);
}

int main() {
	
	class_eq_test<eq1>("eq1");
	class_eq_test<eq2>("eq2");
	
	return 0;
}





