#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
using namespace std ;
#include "Tools.h"
#include "Material.h"
#include "Contact_Law.h"
#include "Node.h"
#include "Gauss.h"
#include "Border.h"
#include "Contact_Element.h"
#include "Body.h"
#include "Spy.h"
#include "IO.h"
#include "Solver.h"
#include "Proximity.h"
#include "Graphic.h"
#include "Monitoring.h"

main(int argc, char **argv)
{
    // CREATE VARIABLES //
    string Simulation_name ;
    int Nb_materials ;
    vector<Material> Materials ;
    int Nb_contact_laws ;
    vector<Contact_law> Contact_laws ;
    string Solver ;
    double Tini, Deltat, Tend, Time ;
    double Target_error, Inv_Target_error, Control_parameter, Accepted_ratio ;
    double Max_mass_scaling, Control_parameter_mass_scaling, Error_factor_mass_scaling, Decrease_factor_mass_scaling ;
    double Save_period, Print_period, Contact_update_period ;
    double Next_save, Next_print, Next_contact_update ;
    int Number_save, Number_print, Number_iteration ;
    double Xmin_period, Xmax_period, Penalty ;
    double Xgravity, Ygravity ;
    double Chains_typical_pressure, Chains_size_ratio ;
    double Fields_xmin, Fields_xmax, Fields_ymin, Fields_ymax, Fields_step, Fields_dist ;
    int Nb_bodies ;
    vector<Body> Bodies ;
    int Nb_monitored ;
    vector<vector<double>> Monitored ;
    int Nb_deactivated ;
    vector<vector<double>> Deactivated ;
    int Nb_spies ;
    vector<Spy> Spies ;
    int Nb_regions = 0 ;
    vector<vector<int>> Regions ;
    vector<int> flags(20) ;
    //int flag_failure = 0 ;
    vector<int> To_Plot(54) ;
    vector<vector<int>> Contacts_Table ;

    // LOAD STATIC DATA //
    Load_static( Simulation_name, Nb_materials, Materials,
                 Nb_contact_laws, Contact_laws, Contacts_Table,
                 Solver, Tini,	Deltat, Tend,
                 Target_error, Inv_Target_error, Control_parameter, Accepted_ratio,
                 Max_mass_scaling, Control_parameter_mass_scaling, Error_factor_mass_scaling, Decrease_factor_mass_scaling,
                 Save_period, Print_period, Contact_update_period,
                 Xmin_period, Xmax_period, Penalty, Xgravity, Ygravity,
                 Chains_typical_pressure, Chains_size_ratio,
                 Fields_xmin, Fields_xmax, Fields_ymin, Fields_ymax, Fields_step, Fields_dist,
                 Nb_monitored, Monitored, Nb_deactivated, Deactivated, Nb_spies, Spies,
                 Nb_regions, Regions, Nb_bodies, Bodies, To_Plot ) ;

    int istart = atoi(argv[1]) ;
    int interval = atoi(argv[2]) ;
    int iend = atoi(argv[3]) ;

    for (int i(istart) ; i<=iend ; i+=interval)
    {
        Number_save = i ;
        for (int j=0 ; j<(int)flags.size() ; j++)
            flags[j] = 0 ;
        flags[7] = 1 ;
        Load_dynamic( Number_save, Number_print,
                      Time, Number_iteration, Deltat, Solver,
                      Next_save, Next_print, Next_contact_update,
                      Xmin_period, Xmax_period,
                      Bodies, Materials, flags ) ;
        if (flags[3]==1)
            Initialize_CZM( Nb_bodies, Bodies, Nb_contact_laws, Contact_laws, flags, Xmin_period, Xmax_period ) ;
        Update_contact_pressures(Nb_bodies, Bodies) ;
        for (int j=0 ; j<Nb_bodies ; j++)
            Bodies[j].Update_current_positions() ;
        Update_proximity(Nb_bodies, Bodies, Xmin_period, Xmax_period, flags) ;
        for (int j=0 ; j<Nb_bodies ; j++)
            Bodies[j].Initialize_contact_forces() ;
        for (int j=0 ; j<Nb_bodies ; j++)
            Bodies[j].Update_contacts(Bodies, Deltat, Xmin_period, Xmax_period) ;
        for (int j=0 ; j<Nb_bodies ; j++)
            Bodies[i].Update_contact_forces(Deltat, Bodies, Nb_contact_laws, Contact_laws, Contacts_Table, Xmin_period, Xmax_period) ;
        for (int j=0 ; j<Nb_bodies ; j++)
            Bodies[j].Send_contact_forces(Bodies) ;
        for (int j=0 ; j<Nb_bodies ; j++)
            Bodies[j].Update_borders(Xmin_period, Xmax_period) ;
        for (int j=0 ; j<Nb_bodies ; j++)
            Bodies[j].Update_damage() ;
        if (flags[11]==1)
            for (int j=0 ; j<Nb_bodies ; j++)
                Bodies[j].Update_initial_damage() ;

        Number_print = i ;
        Write_fields(Nb_bodies, Bodies, Number_iteration, Number_save, Number_print, Time, Xmin_period, Xmax_period, Fields_xmin, Fields_xmax, Fields_ymin, Fields_ymax, Fields_step, Fields_dist) ;
    }
}

