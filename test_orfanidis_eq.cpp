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
unsigned int cfg_data[3] =
{
#include "orfanidis_eq.tstdat"
};

#define READED_SAMPLE_RATE (cfg_data[0])
#define READED_NUMBER_OF_BANDS (cfg_data[1])
#define TEST_VECTORS_NUMBER_OF_SAMPLES (cfg_data[2])

static const double default_max_gain = 1.0;
static const double default_min_gain = 0.0;
static const double default_unit_imp_amp = 1.0;

class test_orfanidis_eq
{
public:
    freq_grid fg;
    eq equalizer;

    test_orfanidis_eq() : equalizer(fg){}
    ~test_orfanidis_eq(){}

    void initialize_data_vectors(vector<vector<eq_single_t> > &vects,
            unsigned int number_of_bands,
            unsigned int data_size) {

        //Cleanup all
        vects.clear();

        //Provide sizes, and set unit impulse
        for(unsigned int i = 0; i < number_of_bands; i++)
        {
            vector<eq_single_t> vect(data_size);
            vect[0] = default_unit_imp_amp;
            vects.push_back(vect);
        }


    }

    bool set_freq_grid(unsigned int number_of_bands)
    {
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
                cout << "Unexpected data read for number of bands : "
                    << READED_NUMBER_OF_BANDS << endl;
                return true;
        }
        return false;
    }

    void process_channel(unsigned int band_number,
            vector<eq_single_t>& test_data)
    {
        //Decrease all gains, without band_number
        vector<eq_single_t> gains_per_band(equalizer.get_number_of_bands(), default_min_gain);
        gains_per_band[band_number] = default_max_gain;

        //Update gains
        equalizer.change_params(gains_per_band);

        //Process data
        for(unsigned int i = 0; i < test_data.size(); i++)
            equalizer.sbs_process(&test_data[i], &test_data[i]);
    }

    void process_eq(vector<vector<eq_single_t> > &test_vectors)
    {
        for(unsigned int i = 0; i < equalizer.get_number_of_bands(); i++)
            process_channel(i, test_vectors[i]);
    }

    void save_cs_files(vector<vector<eq_single_t> > &test_vectors)
    {
        for(unsigned int i = 0; i < equalizer.get_number_of_bands(); i++)
        {
            //Get file name
            std::ostringstream file_name;
            file_name << equalizer.get_string_eq_type();
            file_name << "_" << i << ".tstdat";

            //Write to file
            ofstream data_file;
            data_file.open(file_name.str().c_str());
            if (!data_file.is_open())
                return;
            for(unsigned int j = 0; j < test_vectors[0].size(); j++)
                data_file << test_vectors[i][j] << ", ";
            data_file.close();
        }
    }
};

int main()
{
    //Print input data
	cout << cfg_data[0] << " " << cfg_data[1] << " " << cfg_data[2] << endl;

    test_orfanidis_eq test_eq;

    if(test_eq.set_freq_grid(READED_NUMBER_OF_BANDS)) return 1;

    vector<vector<eq_single_t> > test_vectors;
    test_eq.initialize_data_vectors(test_vectors, READED_NUMBER_OF_BANDS,
                                    TEST_VECTORS_NUMBER_OF_SAMPLES);

	//Butterworth filters based eq calculation, and processing
    test_eq.equalizer.set_eq(test_eq.fg, butterworth);
    test_eq.equalizer.set_sample_rate(READED_SAMPLE_RATE);
    test_eq.process_eq(test_vectors);
    test_eq.save_cs_files(test_vectors);

    test_eq.initialize_data_vectors(test_vectors, READED_NUMBER_OF_BANDS,
                                    TEST_VECTORS_NUMBER_OF_SAMPLES);

    //Chebyshev2
    test_eq.equalizer.set_eq(test_eq.fg, chebyshev1);
    test_eq.equalizer.set_sample_rate(READED_SAMPLE_RATE);
    test_eq.process_eq(test_vectors);
    test_eq.save_cs_files(test_vectors);

    test_eq.initialize_data_vectors(test_vectors, READED_NUMBER_OF_BANDS,
                                     TEST_VECTORS_NUMBER_OF_SAMPLES);

     //Chebyshev2
     test_eq.equalizer.set_eq(test_eq.fg, chebyshev2);
     test_eq.equalizer.set_sample_rate(READED_SAMPLE_RATE);
     test_eq.process_eq(test_vectors);
     test_eq.save_cs_files(test_vectors);

     cout << "Completed" << endl;

	return 0;
}





