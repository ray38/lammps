



void compute(int nneigh, int* nei_list_i, double* nei_list_d)
{
    // function arguments:
    // nneigh: intger, number of neighbors in the neighborlist
    // nei_list_i: in(teger parameters for the neighbor atoms, has the length of "nneigh * 2" assumed the follow format:
    //   [atom1_type, atom1_index, atom2_type, atom2_index, ...]
    //      atom_type: type of atom, could be atomic number, or we can even give subtypes (e.g. C_sp2, C_sp3), as long as we index them correctly
    //      atom_index: index of the neighboring atom (for the right ordering of the resulting derivative matrix)
    // nei_list_d: double parameters for the neighbor atoms, has the length of "nneigh * 4"assumed the following format:
    //   [x_shift_atom1, y_shift_atom1, z_shift_atom1, r_atom1, x_shift_atom2, y_shift_atom2, ....]

    // pre-read / defined
    // ngmp: total number of gmp features per atom
    // square: whether the resulting features should be squared
    // params_i: array of intger parameters for gmp feature definition, 
    //           has length of "ngmp", each element defines the order of the corresponding gmp feature.
    //    [probe1 MCSH order, probe2 MCSH order, ...]
    // parames_d: matrix of double parameters for gmp feature definitions, has the size of "ngmp * 2" follow the following format:
    //    [[probe1 A, probe1 alpha], [probe2 A, probe2 alpha], ...]
    //    * where A, alpha defines the Gaussian probe function A*exp(-alpha * x^2)
    //    * therefore, given the std_dev (sigma) of a gaussian function
    //    * A = (1.0 / (sigma * sqrt(2.0 * PI))) ** 3,
    //    * alpha = 1.0 / (2.0 * sigma * sigma),
    // ngaussian: array of number of gaussian for representing each type of atom, has the length of "num_atom_types"
    //    [num gaussian for H, num gaussians for He, etc]
    // atom_gaussian: matrix to hold the information of primitive gaussian, has the size of "num_atom_types * (max_num_primitive_gaussian * 2)"
    //    [[atom_type1 gaussian1 B, atom_type1 gaussian1 beta, atom_type1 gaussian2 B, atom_type1, gaussian2 beta], ...
    //     [atom_type2 gaussian1 B, atom_type2 gaussian1 beta, atom_type2 gaussian2 B, atom_type2, gaussian2 beta],...]
    //    * B, beta define a specific primitive gaussian, mean the same thing as A, alpha mentioned before, named differently just for easy understanding
    //    * these are read from GMP pseudo potential files.

    // Hard-coded Ta test, small descriptor set
    // 10 sigmas, 3 orders, width 2.0
    int ngmp = 40;
    int square = 1;
    int params_i[40] = {0,0,0,0,0,0,0,0,0,0,
                        1,1,1,1,1,1,1,1,1,1,
                        2,2,2,2,2,2,2,2,2,2,
                        3,3,3,3,3,3,3,3,3,3}; 
    double params_d[40][2] = {{7.936704491780122, 12.5},{0.9920880614725153, 3.125},{0.2939520182140786, 1.3888888888888888},{0.12401100768406441, 0.78125},
                              {0.06349363593424098, 0.5},{0.03674400227675983, 0.3472222222222222},{0.02313908015096246, 0.2551020408163266},
                              {0.015501375960508051, 0.1953125},{0.010887111785706616, 0.15432098765432098},{0.007936704491780123, 0.125},
                              {7.936704491780122, 12.5},{0.9920880614725153, 3.125},{0.2939520182140786, 1.3888888888888888},{0.12401100768406441, 0.78125},
                              {0.06349363593424098, 0.5},{0.03674400227675983, 0.3472222222222222},{0.02313908015096246, 0.2551020408163266},
                              {0.015501375960508051, 0.1953125},{0.010887111785706616, 0.15432098765432098},{0.007936704491780123, 0.125},
                              {7.936704491780122, 12.5},{0.9920880614725153, 3.125},{0.2939520182140786, 1.3888888888888888},{0.12401100768406441, 0.78125},
                              {0.06349363593424098, 0.5},{0.03674400227675983, 0.3472222222222222},{0.02313908015096246, 0.2551020408163266},
                              {0.015501375960508051, 0.1953125},{0.010887111785706616, 0.15432098765432098},{0.007936704491780123, 0.125},
                              {7.936704491780122, 12.5},{0.9920880614725153, 3.125},{0.2939520182140786, 1.3888888888888888},{0.12401100768406441, 0.78125},
                              {0.06349363593424098, 0.5},{0.03674400227675983, 0.3472222222222222},{0.02313908015096246, 0.2551020408163266},
                              {0.015501375960508051, 0.1953125},{0.010887111785706616, 0.15432098765432098},{0.007936704491780123, 0.125},}; 

    int element_index_to_order[] = {0};
    double atom_gaussian[5][2] = {{324.7007374101885,	7.069358713570912},
                                  {0.3192341035226007,	0.5845119391444002},
                                  {2.5692440955326967,	2.183254099956451},
                                  {31.8978274188846,	11.196213264426499},
                                  {-358.43547357714857,	7.539608640456548}};
    int ngaussians[1] = {5};


    // outputs:
    // gmp: array to store the features, has the length of "ngmp"
    // dgmp: matrix for feature derivative, has the size "(ngmp) * (nneigh*3)"

    double* gmp = new double[ngmp];
    double** dgmp = new double[ngmp][nneigh * 3];


    for (int m = 0; m < ngmp; ++m) {
        int mcsh_order = params_i[m];
        int num_groups = get_num_groups(mcsh_order); // defined in gmp_helper
        double A = params_d[m][2], alpha = params_d[m][3];
        double weight = 1.0;
        double sum_square = 0.0;
        for (int group_index = 1; group_index < (num_groups+1); ++group_index){
            SolidGMPFunction gmp_function = get_solid_gmp_function(mcsh_order, group_index); 
            // defined in gmp_sold_harmonics.cpp, get the right gmp function based on MCSH_order and group 
            double group_coefficient = get_group_coefficients(mcsh_order, group_index);
            // defined in gmp_sold_harmonics.cpp, get the corresponding group coefficient
            int mcsh_type = get_mcsh_type(mcsh_order, group_index);
            // defined in gmp_sold_harmonics.cpp, get the type of the mcsh_group (1-membered group, 3-membered group, or 6-membered group)

            if (mcsh_type == 1){

                double sum_miu = 0.0;

                double* sum_dmiu_dxj = new double[nneigh];
                double* sum_dmiu_dyj = new double[nneigh];
                double* sum_dmiu_dzj = new double[nneigh];
                for (int j=0; j<nneigh; j++) {
                    sum_dmiu_dxj[j] = 0.0;
                    sum_dmiu_dyj[j] = 0.0;
                    sum_dmiu_dzj[j] = 0.0;
                }
                double m_desc[1], deriv[3];

                for (int j = 0; j < nneigh; ++j) {

                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        gmp_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                        sum_miu += m_desc[0];
                        sum_dmiu_dxj[j] += deriv[0];
                        sum_dmiu_dyj[j] += deriv[1];
                        sum_dmiu_dzj[j] += deriv[2];
                    }
                }
                sum_square += group_coefficient * sum_miu * sum_miu;

                double dmdx, dmdy, dmdz;
                for (int j = 0; j < nneigh; ++j) {
                    dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                    dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                    dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                    dgmp[m][nei_list_i[j*2 + 1]*3] += dmdx;
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                }

                delete [] sum_dmiu_dxj;
                delete [] sum_dmiu_dyj;
                delete [] sum_dmiu_dzj;
                
            }

            if (mcsh_type == 2){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                double* sum_dmiu1_dxj = new double[nneigh];
                double* sum_dmiu2_dxj = new double[nneigh];
                double* sum_dmiu3_dxj = new double[nneigh];
                double* sum_dmiu1_dyj = new double[nneigh];
                double* sum_dmiu2_dyj = new double[nneigh];
                double* sum_dmiu3_dyj = new double[nneigh];
                double* sum_dmiu1_dzj = new double[nneigh];
                double* sum_dmiu2_dzj = new double[nneigh];
                double* sum_dmiu3_dzj = new double[nneigh];
                for (int j=0; j<nneigh; j++) {
                    sum_dmiu1_dxj[j] = 0.0;
                    sum_dmiu2_dxj[j] = 0.0;
                    sum_dmiu3_dxj[j] = 0.0;
                    sum_dmiu1_dyj[j] = 0.0;
                    sum_dmiu2_dyj[j] = 0.0;
                    sum_dmiu3_dyj[j] = 0.0;
                    sum_dmiu1_dzj[j] = 0.0;
                    sum_dmiu2_dzj[j] = 0.0;
                    sum_dmiu3_dzj[j] = 0.0;
                }

                double miu[3], deriv[9];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        gmp_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                        sum_dmiu1_dxj[j] += deriv[0];
                        sum_dmiu1_dyj[j] += deriv[1];
                        sum_dmiu1_dzj[j] += deriv[2];
                        sum_dmiu2_dxj[j] += deriv[3];
                        sum_dmiu2_dyj[j] += deriv[4];
                        sum_dmiu2_dzj[j] += deriv[5];
                        sum_dmiu3_dxj[j] += deriv[6];
                        sum_dmiu3_dyj[j] += deriv[7];
                        sum_dmiu3_dzj[j] += deriv[8];
                    }
                }
                sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                double dmdx, dmdy, dmdz;
                for (int j = 0; j < nneigh; ++j) {
                    dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                    dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                    dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                    dgmp[m][nei_list_i[j*2 + 1]*3] += dmdx;
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                }

                delete [] sum_dmiu1_dxj;
                delete [] sum_dmiu2_dxj;
                delete [] sum_dmiu3_dxj;
                delete [] sum_dmiu1_dyj;
                delete [] sum_dmiu2_dyj;
                delete [] sum_dmiu3_dyj;
                delete [] sum_dmiu1_dzj;
                delete [] sum_dmiu2_dzj;
                delete [] sum_dmiu3_dzj;
                
            }

            if (mcsh_type == 3){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                double* sum_dmiu1_dxj = new double[nneigh];
                double* sum_dmiu2_dxj = new double[nneigh];
                double* sum_dmiu3_dxj = new double[nneigh];
                double* sum_dmiu4_dxj = new double[nneigh];
                double* sum_dmiu5_dxj = new double[nneigh];
                double* sum_dmiu6_dxj = new double[nneigh];
                double* sum_dmiu1_dyj = new double[nneigh];
                double* sum_dmiu2_dyj = new double[nneigh];
                double* sum_dmiu3_dyj = new double[nneigh];
                double* sum_dmiu4_dyj = new double[nneigh];
                double* sum_dmiu5_dyj = new double[nneigh];
                double* sum_dmiu6_dyj = new double[nneigh];
                double* sum_dmiu1_dzj = new double[nneigh];
                double* sum_dmiu2_dzj = new double[nneigh];
                double* sum_dmiu3_dzj = new double[nneigh];
                double* sum_dmiu4_dzj = new double[nneigh];
                double* sum_dmiu5_dzj = new double[nneigh];
                double* sum_dmiu6_dzj = new double[nneigh];
                for (int j=0; j<nneigh; j++) {
                    sum_dmiu1_dxj[j] = 0.0;
                    sum_dmiu2_dxj[j] = 0.0;
                    sum_dmiu3_dxj[j] = 0.0;
                    sum_dmiu4_dxj[j] = 0.0;
                    sum_dmiu5_dxj[j] = 0.0;
                    sum_dmiu6_dxj[j] = 0.0;
                    sum_dmiu1_dyj[j] = 0.0;
                    sum_dmiu2_dyj[j] = 0.0;
                    sum_dmiu3_dyj[j] = 0.0;
                    sum_dmiu4_dyj[j] = 0.0;
                    sum_dmiu5_dyj[j] = 0.0;
                    sum_dmiu6_dyj[j] = 0.0;
                    sum_dmiu1_dzj[j] = 0.0;
                    sum_dmiu2_dzj[j] = 0.0;
                    sum_dmiu3_dzj[j] = 0.0;
                    sum_dmiu4_dzj[j] = 0.0;
                    sum_dmiu5_dzj[j] = 0.0;
                    sum_dmiu6_dzj[j] = 0.0;
                }

                double miu[6], deriv[18];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        gmp_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                        sum_miu4 += miu[3];
                        sum_miu5 += miu[4];
                        sum_miu6 += miu[5];
                        sum_dmiu1_dxj[j] += deriv[0];
                        sum_dmiu1_dyj[j] += deriv[1];
                        sum_dmiu1_dzj[j] += deriv[2];
                        sum_dmiu2_dxj[j] += deriv[3];
                        sum_dmiu2_dyj[j] += deriv[4];
                        sum_dmiu2_dzj[j] += deriv[5];
                        sum_dmiu3_dxj[j] += deriv[6];
                        sum_dmiu3_dyj[j] += deriv[7];
                        sum_dmiu3_dzj[j] += deriv[8];
                        sum_dmiu4_dxj[j] += deriv[9];
                        sum_dmiu4_dyj[j] += deriv[10];
                        sum_dmiu4_dzj[j] += deriv[11];
                        sum_dmiu5_dxj[j] += deriv[12];
                        sum_dmiu5_dyj[j] += deriv[13];
                        sum_dmiu5_dzj[j] += deriv[14];
                        sum_dmiu6_dxj[j] += deriv[15];
                        sum_dmiu6_dyj[j] += deriv[16];
                        sum_dmiu6_dzj[j] += deriv[17];
                    }
                }
                sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                    sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                double dmdx, dmdy, dmdz;

                for (int j = 0; j < nneigh; ++j) {
                    dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                            sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                            sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                    dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                            sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                            sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                    dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                            sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                            sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                    dgmp[m][nei_list_i[j*2 + 1]*3] += dmdx;
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                }
                
                delete [] sum_dmiu1_dxj;
                delete [] sum_dmiu2_dxj;
                delete [] sum_dmiu3_dxj;
                delete [] sum_dmiu4_dxj;
                delete [] sum_dmiu5_dxj;
                delete [] sum_dmiu6_dxj;
                delete [] sum_dmiu1_dyj;
                delete [] sum_dmiu2_dyj;
                delete [] sum_dmiu3_dyj;
                delete [] sum_dmiu4_dyj;
                delete [] sum_dmiu5_dyj;
                delete [] sum_dmiu6_dyj;
                delete [] sum_dmiu1_dzj;
                delete [] sum_dmiu2_dzj;
                delete [] sum_dmiu3_dzj;
                delete [] sum_dmiu4_dzj;
                delete [] sum_dmiu5_dzj;
                delete [] sum_dmiu6_dzj;
            }
        }
        // sum_square = sum_square * weight;
        if (square != 0){
            gmp[m] = sum_square;
        }
        else {
            double temp = sqrt(sum_square);
            if (fabs(temp) < 1e-8){
                gmp[m] = 0.0;
                for (int j = 0; j < nneigh; ++j) {

                    dgmp[m][nei_list_i[j*2 + 1]*3] = 0.0;
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 1] = 0.0;
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 2] = 0.0;
                }
            }
            else {
                gmp[m] = temp;
                for (int j = 0; j < nneigh; ++j) {

                    dgmp[m][nei_list_i[j*2 + 1]*3] *= (0.5 / temp);
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 1] *= (0.5 / temp);
                    dgmp[m][nei_list_i[j*2 + 1]*3 + 2] *= (0.5 / temp);
                }

            }

        }

    }

}


void read(char* config_filename)
{
    // read the setup files, tentatively designed to have the following structure

    /*
        sigmas: 0.1 0.2 0.3 0.4 0.5
        orders: -1 0 1 2 3
        atom_types: H C O C_sp2
        H: <path>/H_pseduo_density.g
        C: <path>/C_pseduo_density.g
        O: <path>/O_pseduo_density.g
        C_coarse: <path>/C_pseduo_density2.g
    */

    // first two lines defines the setup of the gmp features
    // third line defines available atom types 
    // the rest determines the paths to the GMP pseudo potential files



}



void read_psp(char* psp_filename)
{
    // read a file contains the definition of primitive gaussians for each atom type that has the following structure
    //      B1 beta1
    //      B2 beta2
    //      B3 beta3
    //      ....

    // The number of primitive gaussians are not fixed (could be different for different atom types)
}