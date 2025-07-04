#ifndef DEF_MONITORING
#define DEF_MONITORING

//********************************************//
//** MONITORING ******************************//
//********************************************//

void Monitoring(int Nb_bodies,
                vector<Body>& Bodies,
                double Time,
                double& Deltat,
                int& Number_iteration,
                int& neval,
                double& max_error,
                double& mean_error,
                vector<int>& flags,
                vector<vector<double>>& monitoring,
                int& Nb_monitored,
                vector<vector<double>>& Monitored)
{
    double Dwork = 0. ;
    double Elastic_energy = 0. ;
    double Kinetic_energy = 0. ;
    double px, py, length, vx, vy ;
    vector<double> current_monitoring = { (double)Number_iteration, Time, Deltat, (double)neval, (double)max_error, (double)mean_error } ;

    //cout << "nb monitored " << Nb_monitored << endl ;
    double x ;
    double y ;
    double r ;
    double m ;

    for (int i(0) ; i < Nb_monitored ; i++)
    {
        if (Monitored[i][0] == 0)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
            if (Monitored[i][3] >=0 )
            {
                x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_current ;
                y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_current ;
            }
            else
            {
                if (Bodies[Monitored[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                    {
                        x += Bodies[Monitored[i][2]].nodes[j].x_current * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        y += Bodies[Monitored[i][2]].nodes[j].y_current * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                        m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                    }
                    x /= m ;
                    y /= m ;
                }
                else if (Bodies[Monitored[i][2]].type == "rigid")
                {
                    x = Bodies[Monitored[i][2]].x_current ;
                    y = Bodies[Monitored[i][2]].y_current ;
                    r = Bodies[Monitored[i][2]].r_current ;
                }
            }
            if (Monitored[i][1] == 0)
                current_monitoring.push_back(x) ;
            else if (Monitored[i][1] == 1)
                current_monitoring.push_back(y) ;
            else if (Monitored[i][1] == 2)
                current_monitoring.push_back(r) ;
        }
        else if (Monitored[i][0] == 1)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
            if (Monitored[i][3] >=0 )
            {
                x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_displacement ;
                y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_displacement ;
            }
            else
            {
                if (Bodies[Monitored[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                    {
                        x += Bodies[Monitored[i][2]].nodes[j].x_displacement * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        y += Bodies[Monitored[i][2]].nodes[j].y_displacement * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                        m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                    }
                    x /= m ;
                    y /= m ;
                }
                else if (Bodies[Monitored[i][2]].type == "rigid")
                {
                    x = Bodies[Monitored[i][2]].x_displacement ;
                    y = Bodies[Monitored[i][2]].y_displacement ;
                    r = Bodies[Monitored[i][2]].r_displacement ;
                }
            }
            if (Monitored[i][1] == 0)
                current_monitoring.push_back(x) ;
            else if (Monitored[i][1] == 1)
                current_monitoring.push_back(y) ;
            else if (Monitored[i][1] == 2)
                current_monitoring.push_back(pow(x*x+y*y,0.5)) ;
            else if (Monitored[i][1] == 3)
                current_monitoring.push_back(r) ;
        }
        else if (Monitored[i][0] == 2)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
            if (Monitored[i][3] >=0 )
            {
                x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_velocity ;
                y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_velocity ;
            }
            else
            {
                if (Bodies[Monitored[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                    {
                        x += Bodies[Monitored[i][2]].nodes[j].x_velocity * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        y += Bodies[Monitored[i][2]].nodes[j].y_velocity * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                        m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                    }
                    x /= m ;
                    y /= m ;
                }
                else if (Bodies[Monitored[i][2]].type == "rigid")
                {
                    x = Bodies[Monitored[i][2]].x_velocity ;
                    y = Bodies[Monitored[i][2]].y_velocity ;
                    r = Bodies[Monitored[i][2]].r_velocity ;
                }
            }
            if (Monitored[i][1] == 0)
                current_monitoring.push_back(x) ;
            else if (Monitored[i][1] == 1)
                current_monitoring.push_back(y) ;
            else if (Monitored[i][1] == 2)
                current_monitoring.push_back(pow(x*x+y*y,0.5)) ;
            else if (Monitored[i][1] == 3)
                current_monitoring.push_back(r) ;
        }
        else if (Monitored[i][0] == 3)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
            if (Monitored[i][3] >=0 )
            {
                x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_acceleration ;
                y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_acceleration ;
            }
            else
            {
                if (Bodies[Monitored[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                    {
                        x += Bodies[Monitored[i][2]].nodes[j].x_acceleration * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        y += Bodies[Monitored[i][2]].nodes[j].y_acceleration * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                        m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                    }
                    x /= m ;
                    y /= m ;
                }
                else if (Bodies[Monitored[i][2]].type == "rigid")
                {
                    x = Bodies[Monitored[i][2]].x_acceleration ;
                    y = Bodies[Monitored[i][2]].y_acceleration ;
                    r = Bodies[Monitored[i][2]].r_acceleration ;
                }
            }
            if (Monitored[i][1] == 0)
                current_monitoring.push_back(x) ;
            else if (Monitored[i][1] == 1)
                current_monitoring.push_back(y) ;
            else if (Monitored[i][1] == 2)
                current_monitoring.push_back(pow(x*x+y*y,0.5)) ;
            else if (Monitored[i][1] == 3)
                current_monitoring.push_back(r) ;
        }
        else if (Monitored[i][0] == 4)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << Monitored[i][4] << endl ;
            if (Monitored[i][4] >=0 )
            {
                if (Monitored[i][2] == 0)
                {
                    x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_force ;
                    y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_force ;
                }
                else if (Monitored[i][2] == 1)
                {
                    x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_internal_force ;
                    y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_internal_force ;
                }
                else if (Monitored[i][2] == 2)
                {
                    x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_contact_force ;
                    y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_contact_force ;
                }
                else if (Monitored[i][2] == 3)
                {
                    x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_body_force ;
                    y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_body_force ;
                }
                else if (Monitored[i][2] == 4)
                {
                    x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_dirichlet_force ;
                    y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_dirichlet_force ;
                }
                else if (Monitored[i][2] == 5)
                {
                    x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_neumann_force ;
                    y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_neumann_force ;
                }
                else if (Monitored[i][2] == 6)
                {
                    x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_damping_force ;
                    y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_damping_force ;
                }
                else if (Monitored[i][2] == 7)
                {
                    x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_alid_force ;
                    y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_alid_force ;
                }
            }
            else
            {
                if (Bodies[Monitored[i][3]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                    {
                        if (Monitored[i][2] == 0)
                        {
                            x += Bodies[Monitored[i][3]].nodes[j].x_force ;
                            y += Bodies[Monitored[i][3]].nodes[j].y_force ;
                        }
                        else if (Monitored[i][2] == 1)
                        {
                            x += Bodies[Monitored[i][3]].nodes[j].x_internal_force ;
                            y += Bodies[Monitored[i][3]].nodes[j].y_internal_force ;
                        }
                        else if (Monitored[i][2] == 2)
                        {
                            x += Bodies[Monitored[i][3]].nodes[j].x_contact_force ;
                            y += Bodies[Monitored[i][3]].nodes[j].y_contact_force ;
                        }
                        else if (Monitored[i][2] == 3)
                        {
                            x += Bodies[Monitored[i][3]].nodes[j].x_body_force ;
                            y += Bodies[Monitored[i][3]].nodes[j].y_body_force ;
                        }
                        else if (Monitored[i][2] == 4)
                        {
                            x += Bodies[Monitored[i][3]].nodes[j].x_dirichlet_force ;
                            y += Bodies[Monitored[i][3]].nodes[j].y_dirichlet_force ;
                        }
                        else if (Monitored[i][2] == 5)
                        {
                            x += Bodies[Monitored[i][3]].nodes[j].x_neumann_force ;
                            y += Bodies[Monitored[i][3]].nodes[j].y_neumann_force ;
                        }
                        else if (Monitored[i][2] == 6)
                        {
                            x += Bodies[Monitored[i][3]].nodes[j].x_damping_force ;
                            y += Bodies[Monitored[i][3]].nodes[j].y_damping_force ;
                        }
                        else if (Monitored[i][2] == 7)
                        {
                            x += Bodies[Monitored[i][3]].nodes[j].x_alid_force ;
                            y += Bodies[Monitored[i][3]].nodes[j].y_alid_force ;
                        }
                    }
                }
                else if (Bodies[Monitored[i][3]].type == "rigid")
                {
                    if (Monitored[i][2] == 0)
                    {
                        x = Bodies[Monitored[i][3]].x_force ;
                        y = Bodies[Monitored[i][3]].y_force ;
                        r = Bodies[Monitored[i][3]].r_force ;
                    }
                    else if (Monitored[i][2] == 2)
                    {
                        x = Bodies[Monitored[i][3]].x_contact_force ;
                        y = Bodies[Monitored[i][3]].y_contact_force ;
                        r = Bodies[Monitored[i][3]].r_contact_force ;
                    }
                    else if (Monitored[i][2] == 3)
                    {
                        x = Bodies[Monitored[i][3]].x_body_force ;
                        y = Bodies[Monitored[i][3]].y_body_force ;
                        r = Bodies[Monitored[i][3]].r_body_force ;
                    }
                    else if (Monitored[i][2] == 4)
                    {
                        x = Bodies[Monitored[i][3]].x_dirichlet_force ;
                        y = Bodies[Monitored[i][3]].y_dirichlet_force ;
                        r = Bodies[Monitored[i][3]].r_dirichlet_force ;
                    }
                    else if (Monitored[i][2] == 5)
                    {
                        x = Bodies[Monitored[i][3]].x_neumann_force ;
                        y = Bodies[Monitored[i][3]].y_neumann_force ;
                        r = Bodies[Monitored[i][3]].r_neumann_force ;
                    }
                    else if (Monitored[i][2] == 6)
                    {
                        x = Bodies[Monitored[i][3]].x_damping_force ;
                        y = Bodies[Monitored[i][3]].y_damping_force ;
                        r = Bodies[Monitored[i][3]].r_damping_force ;
                    }
                }
            }
            if (Monitored[i][1] == 0)
                current_monitoring.push_back(x) ;
            else if (Monitored[i][1] == 1)
                current_monitoring.push_back(y) ;
            else if (Monitored[i][1] == 2)
                current_monitoring.push_back(pow(x*x+y*y,0.5)) ;
            else if (Monitored[i][1] == 3)
                current_monitoring.push_back(r) ;
        }
        else if (Monitored[i][0] == 5)
            current_monitoring.push_back(Bodies[Monitored[i][1]].nodes[Monitored[i][2]].jacobian) ;
        else if (Monitored[i][0] == 6)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
            if (Monitored[i][1] == 0)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].Sigmaxx) ;
            else if (Monitored[i][1] == 1)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].Sigmayy) ;
            else if (Monitored[i][1] == 2)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].Sigmaxy) ;
            else if (Monitored[i][1] == 3)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].Sigmazz) ;
            else if (Monitored[i][1] == 4)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaTresca) ;
            else if (Monitored[i][1] == 5)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaVM) ;
            else if (Monitored[i][1] == 6)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaI) ;
            else if (Monitored[i][1] == 7)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaII) ;
            else if (Monitored[i][1] == 8)
                current_monitoring.push_back(Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaIII) ;
        }
        else if (Monitored[i][0] == 7)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << Monitored[i][4] << endl ;
            for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
            {
                if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                {
                    if (Monitored[i][1] == 0)
                        current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].gapn) ;
                    else if (Monitored[i][1] == 1)
                        current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].gapt) ;
                    else if (Monitored[i][1] == 2)
                        current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].xsi) ;
                    else if (Monitored[i][1] == 3)
                        current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].xnorm) ;
                    else if (Monitored[i][1] == 4)
                        current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].ynorm) ;
                    else if (Monitored[i][1] == 5)
                        current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].internal[0]) ;
                    break ;
                }
            }
        }
        else if (Monitored[i][0] == 8)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << endl ;
            Bodies[Monitored[i][1]].Update_damage() ;
            current_monitoring.push_back(Bodies[Monitored[i][1]].damage) ;
        }
        else if (Monitored[i][0] == 9)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << endl ;
            x = 0. ;
            if (Monitored[i][2] >=0 )
            {
                if (Monitored[i][1] == 0 )
                {
                    if (Bodies[Monitored[i][2]].type == "deformable")
                    {
                        for (int j=0 ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                        {
                            x += 0.5 * Bodies[Monitored[i][2]].nodes[j].x_velocity_parameter * Bodies[Monitored[i][2]].nodes[j].x_velocity_parameter * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                            x += 0.5 * Bodies[Monitored[i][2]].nodes[j].y_velocity_parameter * Bodies[Monitored[i][2]].nodes[j].y_velocity_parameter * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                        }
                    }
                    else if (Bodies[Monitored[i][2]].type == "rigid")
                    {
                        x += 0.5 * Bodies[Monitored[i][2]].x_velocity * Bodies[Monitored[i][2]].x_velocity * Bodies[Monitored[i][2]].mass ;
                        x += 0.5 * Bodies[Monitored[i][2]].y_velocity * Bodies[Monitored[i][2]].y_velocity * Bodies[Monitored[i][2]].mass ;
                        x += 0.5 * Bodies[Monitored[i][2]].r_velocity * Bodies[Monitored[i][2]].r_velocity * Bodies[Monitored[i][2]].inertia ;
                    }
                }
                else if (Monitored[i][1] == 1 )
                {
                    if (Bodies[Monitored[i][2]].type == "deformable")
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_regions ; j++)
                        {
                            x += Bodies[Monitored[i][2]].Elastic_energy[j] ;
                        }
                    }
                }
            }
            else if (Monitored[i][2] == -1 )
            {
                if (Monitored[i][1] == 0 )
                {
                    for (int n=0 ; n<Nb_bodies ; n++)
                    {
                        if (Bodies[n].type == "deformable")
                        {
                            for (int j=0 ; j<Bodies[n].nb_nodes ; j++)
                            {
                                x += 0.5 * Bodies[n].nodes[j].x_velocity_parameter * Bodies[n].nodes[j].x_velocity_parameter * Bodies[n].nodes[j].x_mass ;
                                x += 0.5 * Bodies[n].nodes[j].y_velocity_parameter * Bodies[n].nodes[j].y_velocity_parameter * Bodies[n].nodes[j].y_mass ;
                            }
                        }
                        else if (Bodies[n].type == "rigid")
                        {
                            x += 0.5 * Bodies[n].x_velocity * Bodies[n].x_velocity * Bodies[n].mass ;
                            x += 0.5 * Bodies[n].y_velocity * Bodies[n].y_velocity * Bodies[n].mass ;
                            x += 0.5 * Bodies[n].r_velocity * Bodies[n].r_velocity * Bodies[n].inertia ;
                        }
                    }
                }
                else if (Monitored[i][1] == -1 )
                {
                    for (int n=0 ; n<Nb_bodies ; n++)
                    {
                        if (Bodies[n].type == "deformable")
                        {
                            for (int j(0) ; j<Bodies[n].nb_regions ; j++)
                            {
                                x += Bodies[n].Elastic_energy[j] ;
                            }
                        }
                    }
                }
            }
            current_monitoring.push_back(x) ;
        }
        else if (Monitored[i][0] == 10)
        {
            //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << endl ;
            if (Monitored[i][2] >=0 )
            {
                if (Monitored[i][1] == 0 )
                    current_monitoring.push_back(Bodies[Monitored[i][2]].internal_work) ;
                else if (Monitored[i][1] == 1 )
                    current_monitoring.push_back(Bodies[Monitored[i][2]].contact_work) ;
                else if (Monitored[i][1] == 2 )
                    current_monitoring.push_back(Bodies[Monitored[i][2]].body_work) ;
                else if (Monitored[i][1] == 3 )
                    current_monitoring.push_back(Bodies[Monitored[i][2]].dirichlet_work) ;
                else if (Monitored[i][1] == 4 )
                    current_monitoring.push_back(Bodies[Monitored[i][2]].neumann_work) ;
                else if (Monitored[i][1] == 5 )
                    current_monitoring.push_back(Bodies[Monitored[i][2]].damping_work) ;
                else if (Monitored[i][1] == 6 )
                    current_monitoring.push_back(Bodies[Monitored[i][2]].alid_work) ;
            }
            else if (Monitored[i][2] == -1 )
            {
                x = 0. ;
                for (int n=0 ; n<Nb_bodies ; n++)
                {
                    if (Monitored[i][1] == 0 )
                        x += Bodies[n].internal_work ;
                    else if (Monitored[i][1] == 1 )
                        x += Bodies[n].contact_work ;
                    else if (Monitored[i][1] == 2 )
                        x += Bodies[n].body_work ;
                    else if (Monitored[i][1] == 3 )
                        x += Bodies[n].dirichlet_work ;
                    else if (Monitored[i][1] == 4 )
                        x += Bodies[n].neumann_work ;
                    else if (Monitored[i][1] == 5 )
                        x += Bodies[n].damping_work ;
                    else if (Monitored[i][1] == 6 )
                        x += Bodies[n].alid_work ;
                }
                current_monitoring.push_back(x) ;
            }
        }
    }

    if (flags[1]==1)
    {
        for (int i=0 ; i<Nb_bodies ; i++)
        {
            for (int j(0) ; j<Bodies[i].nb_regions ; j++)
            {
                Elastic_energy += Bodies[i].Elastic_energy[j] ;
            }
            for (int j=0 ; j<Bodies[i].nb_borders ; j++)
            {
                for (int k=0 ; k<Bodies[i].borders[j].number_border_nodes ; k++)
                {
                    px = Bodies[i].borders[j].x_bc_pressure[k] ;
                    py = Bodies[i].borders[j].y_bc_pressure[k] ;
                    length = Bodies[i].borders[j].length[k] ;
                    vx = Bodies[i].nodes[Bodies[i].borders[j].border_nodes[k]].x_velocity ;
                    vy = Bodies[i].nodes[Bodies[i].borders[j].border_nodes[k]].y_velocity ;
                    Dwork += ( px * vx + py * vy ) * length * Deltat ;
                }
            }
            for (int j=0 ; j<Bodies[i].nb_nodes ; j++)
            {
                Kinetic_energy += 0.5 * Bodies[i].nodes[j].x_velocity * Bodies[i].nodes[j].x_velocity * Bodies[i].nodes[j].x_mass ;
                Kinetic_energy += 0.5 * Bodies[i].nodes[j].y_velocity * Bodies[i].nodes[j].y_velocity * Bodies[i].nodes[j].y_mass ;
            }
        }
        current_monitoring.push_back(Dwork) ;
        current_monitoring.push_back(Elastic_energy) ;
        current_monitoring.push_back(Kinetic_energy) ;
    }

    if (flags[2]==1)
    {
        double fx0, fy0, fx1, fy1 ;
        for (int i=0 ; i<Bodies[0].borders[0].number_border_nodes ; i++)
        {
            fx0 += Bodies[0].borders[0].x_bc_pressure[i]*Bodies[0].borders[0].length[i] ;
            fy0 += Bodies[0].borders[0].y_bc_pressure[i]*Bodies[0].borders[0].length[i] ;
        }
        for (int i=0 ; i<Bodies[1].borders[1].number_border_nodes ; i++)
        {
            fx1 += Bodies[1].borders[1].x_bc_pressure[i]*Bodies[1].borders[1].length[i] ;
            fy1 += Bodies[1].borders[1].y_bc_pressure[i]*Bodies[1].borders[1].length[i] ;
        }
        current_monitoring.push_back(fx0) ;
        current_monitoring.push_back(fy0) ;
        current_monitoring.push_back(fx1) ;
        current_monitoring.push_back(fy1) ;
    }

    monitoring.push_back(current_monitoring) ;
}



//********************************************//
//** SPYING **********************************//
//********************************************//

void Spying(int Nb_bodies,
            vector<Body>& Bodies,
            double Time,
            double& Deltat,
            int& Number_iteration,
            vector<vector<vector<double>>>& spying,
            int& Nb_spies,
            vector<Spy>& Spies,
            double Xmin_period,
            double Xmax_period,
            int Nb_materials,
            vector<Material> Materials)
{
    bool flag = false ;

    //Recalculate stress on nodes if there is a spy on them
    for (int n(0) ; n < Nb_spies ; n++)
    {
        if ( Time < Spies[n].next_time-1.e-6*Deltat )
            continue ;

        int nb_quantities = Spies[n].nb_quantities ;
        vector<vector<double>> Monitored = Spies[n].quantities ;
        for (int i(0) ; i < nb_quantities ; i++)
        {
            if (Monitored[i][0] == 6)
            {
                #pragma omp parallel
                {
                    #pragma omp parallel
                    for(int i = 0 ; i < Nb_bodies ; i++) //Recalculate stress on each node
                        Bodies[i].Compute_nodal_stresses( Nb_materials, Materials ) ;
                }

                flag = true ; //If we recalculate stress, we do not need to do it again for other spies at this time step
                break ;

            }
            if (flag == true)
                break;
        }
        if (flag == true)
            break;
    }

    //Main loop
    for (int n(0) ; n < Nb_spies ; n++)
    {
        if ( Time < Spies[n].next_time-1.e-6*Deltat )
            continue ;
        Spies[n].next_time += Spies[n].period ;
        int nb_quantities = Spies[n].nb_quantities ;
        vector<vector<double>> Monitored = Spies[n].quantities ;
        vector<double> current_monitoring = { (double)Number_iteration, Time } ;
        double x ;
        double y ;
        double r ;
        double s ;
        double m ;
        double dx ;
        double dy ;
        double mass_x ;
        double mass_y ;
        double mass_tot ;

        for (int i(0) ; i < nb_quantities ; i++)
        {
            if (Monitored[i][0] == 0)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
                if (Monitored[i][1] == 3 )
                {
                    for (int j(0) ; j < Nb_bodies ; j++)
                    {
                        for (int k(0) ; k < Bodies[j].nb_nodes ; k++)
                            current_monitoring.push_back(Bodies[j].nodes[k].x_current) ;
                    }
                    continue ;
                }
                else if (Monitored[i][1] == 4 )
                {
                    for (int j(0) ; j < Nb_bodies ; j++)
                    {
                        for (int k(0) ; k < Bodies[j].nb_nodes ; k++)
                            current_monitoring.push_back(Bodies[j].nodes[k].y_current) ;
                    }
                    continue ;
                }
                if (Monitored[i][3] >=0 )
                {
                    x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_current ;
                    y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_current ;
                }
                else
                {
                    if (Bodies[Monitored[i][2]].type == "deformable")
                    {
                        x = 0. ;
                        y = 0. ;
                        r = 0. ;
                        m = 0. ;
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                        {
                            x += Bodies[Monitored[i][2]].nodes[j].x_current * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                            y += Bodies[Monitored[i][2]].nodes[j].y_current * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                            m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        }
                        x /= m ;
                        y /= m ;
                    }
                    else if (Bodies[Monitored[i][2]].type == "rigid")
                    {
                        x = Bodies[Monitored[i][2]].x_current ;
                        y = Bodies[Monitored[i][2]].y_current ;
                        r = Bodies[Monitored[i][2]].r_current ;
                    }
                }
                if (Monitored[i][1] == 0)
                    current_monitoring.push_back(x) ;
                else if (Monitored[i][1] == 1)
                    current_monitoring.push_back(y) ;
                else if (Monitored[i][1] == 2)
                    current_monitoring.push_back(r) ;
            }
            else if (Monitored[i][0] == 1)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
                if (Monitored[i][3] >=0 )
                {
                    x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_displacement ;
                    y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_displacement ;
                }
                else
                {
                    if (Bodies[Monitored[i][2]].type == "deformable")
                    {
                        x = 0. ;
                        y = 0. ;
                        m = 0. ;
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                        {
                            x += Bodies[Monitored[i][2]].nodes[j].x_displacement * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                            y += Bodies[Monitored[i][2]].nodes[j].y_displacement * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                            m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        }
                        x /= m ;
                        y /= m ;
                    }
                    else if (Bodies[Monitored[i][2]].type == "rigid")
                    {
                        x = Bodies[Monitored[i][2]].x_displacement ;
                        y = Bodies[Monitored[i][2]].y_displacement ;
                        r = Bodies[Monitored[i][2]].r_displacement ;
                    }
                }
                if (Monitored[i][1] == 0)
                    current_monitoring.push_back(x) ;
                else if (Monitored[i][1] == 1)
                    current_monitoring.push_back(y) ;
                else if (Monitored[i][1] == 2)
                    current_monitoring.push_back(pow(x*x+y*y,0.5)) ;
                else if (Monitored[i][1] == 3)
                    current_monitoring.push_back(r) ;
            }
            else if (Monitored[i][0] == 2)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
                if (Monitored[i][3] >=0 )
                {
                    x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_velocity ;
                    y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_velocity ;
                }
                else
                {
                    if (Bodies[Monitored[i][2]].type == "deformable")
                    {
                        x = 0. ;
                        y = 0. ;
                        m = 0. ;
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                        {
                            x += Bodies[Monitored[i][2]].nodes[j].x_velocity * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                            y += Bodies[Monitored[i][2]].nodes[j].y_velocity * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                            m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        }
                        x /= m ;
                        y /= m ;
                    }
                    else if (Bodies[Monitored[i][2]].type == "rigid")
                    {
                        x = Bodies[Monitored[i][2]].x_velocity ;
                        y = Bodies[Monitored[i][2]].y_velocity ;
                        r = Bodies[Monitored[i][2]].r_velocity ;
                    }
                }
                if (Monitored[i][1] == 0)
                    current_monitoring.push_back(x) ;
                else if (Monitored[i][1] == 1)
                    current_monitoring.push_back(y) ;
                else if (Monitored[i][1] == 2)
                    current_monitoring.push_back(pow(x*x+y*y,0.5)) ;
                else if (Monitored[i][1] == 3)
                    current_monitoring.push_back(r) ;
            }
            else if (Monitored[i][0] == 3)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
                if (Monitored[i][3] >=0 )
                {
                    x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_acceleration ;
                    y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_acceleration ;
                }
                else
                {
                    if (Bodies[Monitored[i][2]].type == "deformable")
                    {
                        x = 0. ;
                        y = 0. ;
                        m = 0. ;
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                        {
                            x += Bodies[Monitored[i][2]].nodes[j].x_acceleration * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                            y += Bodies[Monitored[i][2]].nodes[j].y_acceleration * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                            m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        }
                        x /= m ;
                        y /= m ;
                    }
                    else if (Bodies[Monitored[i][2]].type == "rigid")
                    {
                        x = Bodies[Monitored[i][2]].x_acceleration ;
                        y = Bodies[Monitored[i][2]].y_acceleration ;
                        r = Bodies[Monitored[i][2]].r_acceleration ;
                    }
                }
                if (Monitored[i][1] == 0)
                    current_monitoring.push_back(x) ;
                else if (Monitored[i][1] == 1)
                    current_monitoring.push_back(y) ;
                else if (Monitored[i][1] == 2)
                    current_monitoring.push_back(pow(x*x+y*y,0.5)) ;
                else if (Monitored[i][1] == 3)
                    current_monitoring.push_back(r) ;
            }
            else if (Monitored[i][0] == 4)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << Monitored[i][4] << endl ;
                if (Monitored[i][4] >=0 )
                {
                    if (Monitored[i][2] == 0)
                    {
                        x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_force ;
                        y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_force ;
                    }
                    else if (Monitored[i][2] == 1)
                    {
                        x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_internal_force ;
                        y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_internal_force ;
                    }
                    else if (Monitored[i][2] == 2)
                    {
                        x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_contact_force ;
                        y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_contact_force ;
                    }
                    else if (Monitored[i][2] == 3)
                    {
                        x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_body_force ;
                        y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_body_force ;
                    }
                    else if (Monitored[i][2] == 4)
                    {
                        x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_dirichlet_force ;
                        y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_dirichlet_force ;
                    }
                    else if (Monitored[i][2] == 5)
                    {
                        x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_neumann_force ;
                        y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_neumann_force ;
                    }
                    else if (Monitored[i][2] == 6)
                    {
                        x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_damping_force ;
                        y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_damping_force ;
                    }
                    else if (Monitored[i][2] == 7)
                    {
                        x = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].x_alid_force ;
                        y = Bodies[Monitored[i][3]].nodes[Monitored[i][4]].y_alid_force ;
                    }
                }
                else
                {
                    if (Bodies[Monitored[i][3]].type == "deformable")
                    {
                        x = 0. ;
                        y = 0. ;
                        for (int j(0) ; j<Bodies[Monitored[i][3]].nb_nodes ; j++)
                        {
                            if (Monitored[i][2] == 0)
                            {
                                x += Bodies[Monitored[i][3]].nodes[j].x_force ;
                                y += Bodies[Monitored[i][3]].nodes[j].y_force ;
                            }
                            else if (Monitored[i][2] == 1)
                            {
                                x += Bodies[Monitored[i][3]].nodes[j].x_internal_force ;
                                y += Bodies[Monitored[i][3]].nodes[j].y_internal_force ;
                            }
                            else if (Monitored[i][2] == 2)
                            {
                                x += Bodies[Monitored[i][3]].nodes[j].x_contact_force ;
                                y += Bodies[Monitored[i][3]].nodes[j].y_contact_force ;
                            }
                            else if (Monitored[i][2] == 3)
                            {
                                x += Bodies[Monitored[i][3]].nodes[j].x_body_force ;
                                y += Bodies[Monitored[i][3]].nodes[j].y_body_force ;
                            }
                            else if (Monitored[i][2] == 4)
                            {
                                x += Bodies[Monitored[i][3]].nodes[j].x_dirichlet_force ;
                                y += Bodies[Monitored[i][3]].nodes[j].y_dirichlet_force ;
                            }
                            else if (Monitored[i][2] == 5)
                            {
                                x += Bodies[Monitored[i][3]].nodes[j].x_neumann_force ;
                                y += Bodies[Monitored[i][3]].nodes[j].y_neumann_force ;
                            }
                            else if (Monitored[i][2] == 6)
                            {
                                x += Bodies[Monitored[i][3]].nodes[j].x_damping_force ;
                                y += Bodies[Monitored[i][3]].nodes[j].y_damping_force ;
                            }
                            else if (Monitored[i][2] == 7)
                            {
                                x += Bodies[Monitored[i][3]].nodes[j].x_alid_force ;
                                y += Bodies[Monitored[i][3]].nodes[j].y_alid_force ;
                            }
                        }
                    }
                    else if (Bodies[Monitored[i][3]].type == "rigid")
                    {
                        if (Monitored[i][2] == 0)
                        {
                            x = Bodies[Monitored[i][3]].x_force ;
                            y = Bodies[Monitored[i][3]].y_force ;
                            r = Bodies[Monitored[i][3]].r_force ;
                        }
                        else if (Monitored[i][2] == 2)
                        {
                            x = Bodies[Monitored[i][3]].x_contact_force ;
                            y = Bodies[Monitored[i][3]].y_contact_force ;
                            r = Bodies[Monitored[i][3]].r_contact_force ;
                        }
                        else if (Monitored[i][2] == 3)
                        {
                            x = Bodies[Monitored[i][3]].x_body_force ;
                            y = Bodies[Monitored[i][3]].y_body_force ;
                            r = Bodies[Monitored[i][3]].r_body_force ;
                        }
                        else if (Monitored[i][2] == 4)
                        {
                            x = Bodies[Monitored[i][3]].x_dirichlet_force ;
                            y = Bodies[Monitored[i][3]].y_dirichlet_force ;
                            r = Bodies[Monitored[i][3]].r_dirichlet_force ;
                        }
                        else if (Monitored[i][2] == 5)
                        {
                            x = Bodies[Monitored[i][3]].x_neumann_force ;
                            y = Bodies[Monitored[i][3]].y_neumann_force ;
                            r = Bodies[Monitored[i][3]].r_neumann_force ;
                        }
                        else if (Monitored[i][2] == 6)
                        {
                            x = Bodies[Monitored[i][3]].x_damping_force ;
                            y = Bodies[Monitored[i][3]].y_damping_force ;
                            r = Bodies[Monitored[i][3]].r_damping_force ;
                        }
                    }
                }
                if (Monitored[i][1] == 0)
                    current_monitoring.push_back(x) ;
                else if (Monitored[i][1] == 1)
                    current_monitoring.push_back(y) ;
                else if (Monitored[i][1] == 2)
                    current_monitoring.push_back(pow(x*x+y*y,0.5)) ;
                else if (Monitored[i][1] == 3)
                    current_monitoring.push_back(r) ;
            }
            else if (Monitored[i][0] == 5)
                current_monitoring.push_back(Bodies[Monitored[i][1]].nodes[Monitored[i][2]].jacobian) ;

            else if (Monitored[i][0] == 6)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
                //Bodies[i].Compute_nodal_stresses( Nb_materials, Materials ) ;

                double s, m ;
                if  (Monitored[i][3] >= 0)
                {
                    if (Monitored[i][1] == 0)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].Sigmaxx ;
                    else if (Monitored[i][1] == 1)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].Sigmayy ;
                    else if (Monitored[i][1] == 2)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].Sigmaxy ;
                    else if (Monitored[i][1] == 3)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].Sigmazz ;
                    else if (Monitored[i][1] == 4)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaTresca ;
                    else if (Monitored[i][1] == 5)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaVM ;
                    else if (Monitored[i][1] == 6)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaI ;
                    else if (Monitored[i][1] == 7)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaII ;
                    else if (Monitored[i][1] == 8)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaIII ;
                    else if (Monitored[i][1] == 9)
                        s = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].SigmaSph ;
                }
                else
                {
                    s = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                    {
                        if (Monitored[i][1] == 0)
                            s += Bodies[Monitored[i][2]].nodes[j].Sigmaxx * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 1)
                            s += Bodies[Monitored[i][2]].nodes[j].Sigmayy * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 2)
                            s += Bodies[Monitored[i][2]].nodes[j].Sigmaxy * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 3)
                            s += Bodies[Monitored[i][2]].nodes[j].Sigmazz * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 4)
                            s += Bodies[Monitored[i][2]].nodes[j].SigmaTresca * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 5)
                            s += Bodies[Monitored[i][2]].nodes[j].SigmaVM * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 6)
                            s += Bodies[Monitored[i][2]].nodes[j].SigmaI * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 7)
                            s += Bodies[Monitored[i][2]].nodes[j].SigmaII * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 8)
                            s += Bodies[Monitored[i][2]].nodes[j].SigmaIII * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        else if (Monitored[i][1] == 9)
                            s += Bodies[Monitored[i][2]].nodes[j].SigmaSph * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                        m += Bodies[Monitored[i][2]].nodes[j].x_mass ;
                    }
                    s = s / m ;
                }
                current_monitoring.push_back(s) ;
            }
            else if (Monitored[i][0] == 7)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << Monitored[i][4] << endl ;
                if (Monitored[i][1] == 0)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].gapn ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].gapn ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 1)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].gapt ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].gapt ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 2)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].xsi ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].xsi ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 3)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].xnorm ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].xnorm ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 4)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].ynorm ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].ynorm ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 5)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].internal[0] ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].internal[0] ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 6)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].fx ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].fx ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 7)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].fy ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].fy ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 8)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].length ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].length ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 9)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].contact_elements[j].bodyM ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].contact_elements[j].bodyM ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 10)
                    current_monitoring.push_back(Bodies[Monitored[i][2]].borders[Monitored[i][3]].x_contact_pressure[Monitored[i][4]]) ;
                else if (Monitored[i][1] == 11)
                    current_monitoring.push_back(Bodies[Monitored[i][2]].borders[Monitored[i][3]].y_contact_pressure[Monitored[i][4]]) ;
                else if (Monitored[i][1] == 12)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Monitored[i][2] ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Monitored[i][2] ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 13)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].nodes[Bodies[Monitored[i][2]].contact_elements[j].nodeS].x_current ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].nodes[Bodies[Monitored[i][2]].contact_elements[j].nodeS].x_current ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 14)
                {
                    double q = 0. ;
                    if ((Monitored[i][3]>=0) & (Monitored[i][4]>=0))
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                            {
                                q = Bodies[Monitored[i][2]].nodes[Bodies[Monitored[i][2]].contact_elements[j].nodeS].y_current ;
                                break ;
                            }
                        }
                        current_monitoring.push_back(q) ;
                    }
                    else
                    {
                        for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                        {
                            q = Bodies[Monitored[i][2]].nodes[Bodies[Monitored[i][2]].contact_elements[j].nodeS].y_current ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                /*
                int flagexist(0) ;
                for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                {
                    if (Bodies[Monitored[i][2]].contact_elements[j].borderS == Monitored[i][3] && Bodies[Monitored[i][2]].contact_elements[j].border_nodeS == Monitored[i][4])
                    {
                        if (Monitored[i][1] == 0)           current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].gapn) ;
                        else if (Monitored[i][1] == 1)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].gapt) ;
                        else if (Monitored[i][1] == 2)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].xsi) ;
                        else if (Monitored[i][1] == 3)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].xnorm) ;
                        else if (Monitored[i][1] == 4)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].ynorm) ;
                        else if (Monitored[i][1] == 5)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].internal[0]) ;
                        else if (Monitored[i][1] == 6)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].fx) ;
                        else if (Monitored[i][1] == 7)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].fy) ;
                        else if (Monitored[i][1] == 8)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].length) ;
                        else if (Monitored[i][1] == 9)      current_monitoring.push_back(Bodies[Monitored[i][2]].contact_elements[j].bodyM) ;
                        flagexist = 1 ;
                        break ;
                    }
                }
                if (flagexist == 0) current_monitoring.push_back(0.) ;
                */
            }
            else if (Monitored[i][0] == 8)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << endl ;
                //Bodies[Monitored[i][1]].Update_damage() ;
                if (Monitored[i][1] == 0 )       current_monitoring.push_back(Bodies[Monitored[i][2]].initial_damage) ;
                else if (Monitored[i][1] == 1 )  current_monitoring.push_back(Bodies[Monitored[i][2]].damage) ;
                else if (Monitored[i][1] == 2 )  current_monitoring.push_back((Bodies[Monitored[i][2]].damage - Bodies[Monitored[i][2]].initial_damage) / (1. - Bodies[Monitored[i][2]].initial_damage)) ;
            }
            else if (Monitored[i][0] == 9)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << endl ;
                x = 0. ;
                if (Monitored[i][2] >=0 )
                {
                    if (Monitored[i][1] == 0 )
                    {
                        if (Bodies[Monitored[i][2]].type == "deformable")
                        {
                            for (int j=0 ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                            {
                                x += 0.5 * Bodies[Monitored[i][2]].nodes[j].x_velocity_parameter * Bodies[Monitored[i][2]].nodes[j].x_velocity_parameter * Bodies[Monitored[i][2]].nodes[j].x_mass ;
                                x += 0.5 * Bodies[Monitored[i][2]].nodes[j].y_velocity_parameter * Bodies[Monitored[i][2]].nodes[j].y_velocity_parameter * Bodies[Monitored[i][2]].nodes[j].y_mass ;
                            }
                        }
                        else if (Bodies[Monitored[i][2]].type == "rigid")
                        {
                            x += 0.5 * Bodies[Monitored[i][2]].x_velocity * Bodies[Monitored[i][2]].x_velocity * Bodies[Monitored[i][2]].mass ;
                            x += 0.5 * Bodies[Monitored[i][2]].y_velocity * Bodies[Monitored[i][2]].y_velocity * Bodies[Monitored[i][2]].mass ;
                            x += 0.5 * Bodies[Monitored[i][2]].r_velocity * Bodies[Monitored[i][2]].r_velocity * Bodies[Monitored[i][2]].inertia ;
                        }
                    }
                    else if (Monitored[i][1] == -1 )
                    {
                        if (Bodies[Monitored[i][2]].type == "deformable")
                        {
                            for (int j(0) ; j<Bodies[Monitored[i][2]].nb_regions ; j++)
                            {
                                x += Bodies[Monitored[i][2]].Elastic_energy[j] ;
                            }
                        }
                    }
                }
                else if (Monitored[i][2] == -1 )
                {
                    if (Monitored[i][1] == 0 )
                    {
                        for (int n=0 ; n<Nb_bodies ; n++)
                        {
                            if (Bodies[n].type == "deformable")
                            {
                                for (int j=0 ; j<Bodies[n].nb_nodes ; j++)
                                {
                                    x += 0.5 * Bodies[n].nodes[j].x_velocity_parameter * Bodies[n].nodes[j].x_velocity_parameter * Bodies[n].nodes[j].x_mass ;
                                    x += 0.5 * Bodies[n].nodes[j].y_velocity_parameter * Bodies[n].nodes[j].y_velocity_parameter * Bodies[n].nodes[j].y_mass ;
                                }
                            }
                            else if (Bodies[n].type == "rigid")
                            {
                                x += 0.5 * Bodies[n].x_velocity * Bodies[n].x_velocity * Bodies[n].mass ;
                                x += 0.5 * Bodies[n].y_velocity * Bodies[n].y_velocity * Bodies[n].mass ;
                                x += 0.5 * Bodies[n].r_velocity * Bodies[n].r_velocity * Bodies[n].inertia ;
                            }
                        }
                    }
                    else if (Monitored[i][1] == -1 )
                    {
                        for (int n=0 ; n<Nb_bodies ; n++)
                        {
                            if (Bodies[n].type == "deformable")
                            {
                                for (int j(0) ; j<Bodies[n].nb_regions ; j++)
                                {
                                    x += Bodies[n].Elastic_energy[j] ;
                                }
                            }
                        }
                    }
                }
                current_monitoring.push_back(x) ;
            }
            else if (Monitored[i][0] == 10)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << endl ;
                if (Monitored[i][2] >=0 )
                {
                    if (Monitored[i][1] == 0 )
                        current_monitoring.push_back(Bodies[Monitored[i][2]].internal_work) ;
                    else if (Monitored[i][1] == 1 )
                        current_monitoring.push_back(Bodies[Monitored[i][2]].contact_work) ;
                    else if (Monitored[i][1] == 2 )
                        current_monitoring.push_back(Bodies[Monitored[i][2]].body_work) ;
                    else if (Monitored[i][1] == 3 )
                        current_monitoring.push_back(Bodies[Monitored[i][2]].dirichlet_work) ;
                    else if (Monitored[i][1] == 4 )
                        current_monitoring.push_back(Bodies[Monitored[i][2]].neumann_work) ;
                    else if (Monitored[i][1] == 5 )
                        current_monitoring.push_back(Bodies[Monitored[i][2]].damping_work) ;
                    else if (Monitored[i][1] == 6 )
                        current_monitoring.push_back(Bodies[Monitored[i][2]].alid_work) ;
                }
                else if (Monitored[i][2] == -1 )
                {
                    x = 0. ;
                    for (int n=0 ; n<Nb_bodies ; n++)
                    {
                        if (Monitored[i][1] == 0 )
                            x += Bodies[n].internal_work ;
                        else if (Monitored[i][1] == 1 )
                            x += Bodies[n].contact_work ;
                        else if (Monitored[i][1] == 2 )
                            x += Bodies[n].body_work ;
                        else if (Monitored[i][1] == 3 )
                            x += Bodies[n].dirichlet_work ;
                        else if (Monitored[i][1] == 4 )
                            x += Bodies[n].neumann_work ;
                        else if (Monitored[i][1] == 5 )
                            x += Bodies[n].damping_work ;
                        else if (Monitored[i][1] == 6 )
                            x += Bodies[n].alid_work ;
                    }
                    current_monitoring.push_back(x) ;
                }
            }
            else if (Monitored[i][0] == 11) //Spy ERROR
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << endl ;
                if (Monitored[i][3] >=0 )
                {
                    x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_error ; //Error along x for node i
                    y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_error ;
                    s = pow(x*x + y*y, 0.5) ; //Norm error along x and y for node i
                    r = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].error_norm ; //Error normalized
                    m = r ; //useless, there is no "max" on a single node, this line if for avoiding a bug
                }
                else if (Monitored[i][3] == -1 )
                {
                    x = 0. ;
                    y = 0. ;
                    s = 0. ;
                    r = 0. ;
                    m = 0. ;
                    for (int n = 0 ; n<Bodies[Monitored[i][2]].nb_nodes ; n++)
                    {
                        if (Monitored[i][1] == 0 )
                            x += Bodies[Monitored[i][2]].nodes[n].x_error ;
                        else if (Monitored[i][1] == 1 )
                            y += Bodies[Monitored[i][2]].nodes[n].y_error ;
                        else if (Monitored[i][1] == 2 )
                            s += pow(Bodies[Monitored[i][2]].nodes[n].x_error*Bodies[Monitored[i][2]].nodes[n].x_error +
                                     Bodies[Monitored[i][2]].nodes[n].y_error*Bodies[Monitored[i][2]].nodes[n].y_error, 0.5) ;
                        else if (Monitored[i][1] == 3 )
                            r += Bodies[Monitored[i][2]].nodes[n].error_norm ;
                        else if (Monitored[i][1] == 4 )
                        {
                            if (Bodies[Monitored[i][2]].nodes[n].error_norm > m)
                                m = Bodies[Monitored[i][2]].nodes[n].error_norm ;
                        }
                    }
                    x /= Bodies[Monitored[i][2]].nb_nodes ; //These quantities are averaged on all nodes, except m for max
                    y /= Bodies[Monitored[i][2]].nb_nodes ;
                    s /= Bodies[Monitored[i][2]].nb_nodes ;
                    r /= Bodies[Monitored[i][2]].nb_nodes ;
                }

                if (Monitored[i][1] == 0)
                    current_monitoring.push_back(x) ;
                else if (Monitored[i][1] == 1)
                    current_monitoring.push_back(y) ;
                else if (Monitored[i][1] == 2)
                    current_monitoring.push_back(s) ;
                else if (Monitored[i][1] == 3)
                    current_monitoring.push_back(r) ;
                else if (Monitored[i][1] == 4)
                    current_monitoring.push_back(m) ;
            }
            else if (Monitored[i][0] == 12) // Spy MASS SCALING
            {
                if (Monitored[i][3] >=0 ) //On a node
                {
                    x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_factor_mass_scaling ;
                    y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_factor_mass_scaling ;
                    r = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].factor_mass_scaling ;

                    dx = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].delta_x_factor_mass_scaling ;
                    dy = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].delta_y_factor_mass_scaling ;

                    mass_x = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_factor_mass_scaling * Bodies[Monitored[i][2]].nodes[Monitored[i][3]].x_mass ;
                    mass_y = Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_factor_mass_scaling * Bodies[Monitored[i][2]].nodes[Monitored[i][3]].y_mass ;
                    mass_tot = sqrt(mass_x * mass_x + mass_y * mass_y) ;
                }
                else if (Monitored[i][3] == -1 ) //On every node of one body
                {
                    x = 0. ;
                    y = 0. ;
                    r = 0. ;
                    dx = 0. ;
                    dy = 0. ;
                    mass_x = 0. ;
                    mass_y = 0. ;
                    mass_tot = 0. ;

                    for (int n=0 ; n<Bodies[Monitored[i][2]].nb_nodes ; n++) //Main loop
                    {
                        if (Monitored[i][1] == 0 )
                            x += Bodies[Monitored[i][2]].nodes[n].x_factor_mass_scaling ;
                        else if (Monitored[i][1] == 1 )
                            y += Bodies[Monitored[i][2]].nodes[n].y_factor_mass_scaling ;
                        else if (Monitored[i][1] == 2 )
                            r += Bodies[Monitored[i][2]].nodes[n].factor_mass_scaling ;

                        else if (Monitored[i][1] == 3 )
                            dx += Bodies[Monitored[i][2]].nodes[n].delta_x_factor_mass_scaling ;
                        else if (Monitored[i][1] == 4)
                            dy += Bodies[Monitored[i][2]].nodes[n].delta_y_factor_mass_scaling ;
                        else if (Monitored[i][1] == 5)
                        {
                            dx += Bodies[Monitored[i][2]].nodes[n].delta_x_factor_mass_scaling ;
                            dy += Bodies[Monitored[i][2]].nodes[n].delta_y_factor_mass_scaling ;
                        }

                        else if (Monitored[i][1] == 6)
                            mass_x += Bodies[Monitored[i][2]].nodes[n].x_factor_mass_scaling * Bodies[Monitored[i][2]].nodes[n].x_mass ;
                        else if (Monitored[i][1] == 7)
                            mass_y += Bodies[Monitored[i][2]].nodes[n].y_factor_mass_scaling * Bodies[Monitored[i][2]].nodes[n].y_mass ;
                        else if (Monitored[i][1] == 8)
                        {
                            mass_x = Bodies[Monitored[i][2]].nodes[n].x_factor_mass_scaling * Bodies[Monitored[i][2]].nodes[n].x_mass ;
                            mass_y = Bodies[Monitored[i][2]].nodes[n].y_factor_mass_scaling * Bodies[Monitored[i][2]].nodes[n].y_mass ;
                            mass_tot += sqrt(mass_x * mass_x + mass_y * mass_y) ;
                        }
                    }
                }

                if (Monitored[i][1] == 0)
                    current_monitoring.push_back(x) ;
                else if (Monitored[i][1] == 1)
                    current_monitoring.push_back(y) ;
                else if (Monitored[i][1] == 2)
                    current_monitoring.push_back(r) ;
                else if (Monitored[i][1] == 3)
                    current_monitoring.push_back(dx) ;
                else if (Monitored[i][1] == 4)
                    current_monitoring.push_back(dy) ;
                else if (Monitored[i][1] == 5)
                    current_monitoring.push_back((dx + dy) * 0.5) ;
                else if (Monitored[i][1] == 6)
                    current_monitoring.push_back(mass_x) ;
                else if (Monitored[i][1] == 7)
                    current_monitoring.push_back(mass_y) ;
                else if (Monitored[i][1] == 8)
                    current_monitoring.push_back(mass_tot) ;
            }
            else if (Monitored[i][0] == 13)
            {
                int n = Monitored[i][2] ;
                double interval = 3.141592653589793 / n ;
                vector<double> Bins(n) ;
                int inbin ;
                double period = Xmax_period - Xmin_period ;
                double angle = 0. ;
                double xs, ys ;
                for (int k=0 ; k<n ; k++)
                    Bins[k] = 0. ;
                for (int k=0 ; k<Nb_bodies ; k++)
                {
                    if (Bodies[k].status == "inactive")
                        continue ;
                    xs = Bodies[k].x_current - period * floor( ( Bodies[k].x_current - Xmin_period ) / period ) ;
                    ys = Bodies[k].y_current ;
                    if ((xs<Monitored[i][3]) || (xs>Monitored[i][4]) || (ys<Monitored[i][5]) || (ys>Monitored[i][6]))
                        continue ;
                    for (int j=0 ; j<Bodies[k].nb_contact_elements ; j++)
                    {
                        if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                            continue ;
                        if ( (Monitored[i][1] == 2) || (Monitored[i][1] == 3) )
                        {
                            if (Bodies[k].contact_elements[j].fx * Bodies[k].contact_elements[j].xnorm + Bodies[k].contact_elements[j].fy * Bodies[k].contact_elements[j].ynorm < 0.)
                                continue ;
                        }
                        if ( (Monitored[i][1] == 4) || (Monitored[i][1] == 5) )
                        {
                            if (Bodies[k].contact_elements[j].fx * Bodies[k].contact_elements[j].xnorm + Bodies[k].contact_elements[j].fy * Bodies[k].contact_elements[j].ynorm > 0.)
                                continue ;
                        }
                        if ( (Monitored[i][1] == 0) || (Monitored[i][1] == 2) || (Monitored[i][1] == 4) )
                            angle = atan(Bodies[k].contact_elements[j].ynorm/Bodies[k].contact_elements[j].xnorm) ;
                        else if ( (Monitored[i][1] == 1) || (Monitored[i][1] == 3) || (Monitored[i][1] == 5) )
                            angle = atan(Bodies[k].contact_elements[j].fy/Bodies[k].contact_elements[j].fx) ;
                        if (angle<0.)
                            angle = angle + 3.141592653589793 ;
                        inbin = (int) floor(angle / interval) ;
                        Bins[inbin] = Bins[inbin] + 1. ;
                    }
                }
                for (int k=0 ; k<n ; k++)
                    current_monitoring.push_back(Bins[k]) ;
            }
            else if (Monitored[i][0] == 14)
            {
                if (Monitored[i][1] == 0)
                {
                    if (Bodies[Monitored[i][2]].type=="inactive")
                        current_monitoring.push_back(0) ;
                    else if (Bodies[Monitored[i][2]].type=="active")
                        current_monitoring.push_back(1) ;
                }
                if (Monitored[i][1] == 1)
                {
                    if (Bodies[Monitored[i][2]].type=="rigid")
                        current_monitoring.push_back(0) ;
                    else if (Bodies[Monitored[i][2]].type=="deformable")
                        current_monitoring.push_back(1) ;
                }
                if (Monitored[i][1] == 2)
                {
                    current_monitoring.push_back(Bodies[Monitored[i][2]].nb_active_contacts) ;
                }
                if (Monitored[i][1] == 3)
                {
                    current_monitoring.push_back(Bodies[Monitored[i][2]].nb_contacting_bodies) ;
                }
                else if (Monitored[i][1] == 4)
                {
                    double l = 0. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_borders ; j++)
                    {
                        for (int k(0) ; k<Bodies[Monitored[i][2]].borders[j].number_border_nodes ; k++)
                        {
                            l += Bodies[Monitored[i][2]].borders[j].length[k] ;
                        }
                    }
                    current_monitoring.push_back(l) ;
                }
                else if (Monitored[i][1] == 5)
                {
                    double l = 0. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_contact_elements ; j++)
                    {
                        if (Bodies[Monitored[i][2]].contact_elements[j].fx != 0. && Bodies[Monitored[i][2]].contact_elements[j].fy != 0.)
                        {
                            l += Bodies[Monitored[i][2]].contact_elements[j].length ;
                        }
                    }
                    current_monitoring.push_back(l) ;
                }
                else if (Monitored[i][1] == 6)
                {
                    double a = 0. ;
                    double x1, x2, y1, y2 ;
                    int n1, n2 ;
                    //double y0 = 1.e10 ;
                    //for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                    //{
                    //    if (Bodies[Monitored[i][2]].nodes[j].y_current < y0) y0 = Bodies[Monitored[i][2]].nodes[j].y_current ;
                    //}
                    //y0 -= 1. ;
                    for (int j(0) ; j<Bodies[Monitored[i][2]].nb_borders ; j++)
                    {
                        for (int k(0) ; k<Bodies[Monitored[i][2]].borders[j].number_border_nodes - 1 ; k++)
                        {
                            n1 = Bodies[Monitored[i][2]].borders[j].border_nodes[k] ;
                            n2 = Bodies[Monitored[i][2]].borders[j].border_nodes[k + 1] ;
                            x1 = Bodies[Monitored[i][2]].nodes[n1].x_current ;
                            y1 = Bodies[Monitored[i][2]].nodes[n1].y_current ;// - y0 ;
                            x2 = Bodies[Monitored[i][2]].nodes[n2].x_current ;
                            y2 = Bodies[Monitored[i][2]].nodes[n2].y_current ;// - y0 ;
                            a += 0.5 * ( y1 + y2 ) * ( x1 - x2 ) ;
                        }
                    }
                    n1 = Bodies[Monitored[i][2]].borders[Bodies[Monitored[i][2]].nb_borders-1].border_nodes[Bodies[Monitored[i][2]].borders[Bodies[Monitored[i][2]].nb_borders-1].number_border_nodes-1] ;
                    n2 = Bodies[Monitored[i][2]].borders[0].border_nodes[0] ;
                    x1 = Bodies[Monitored[i][2]].nodes[n1].x_current ;
                    y1 = Bodies[Monitored[i][2]].nodes[n1].y_current ;// - y0 ;
                    x2 = Bodies[Monitored[i][2]].nodes[n2].x_current ;
                    y2 = Bodies[Monitored[i][2]].nodes[n2].y_current ;// - y0 ;
                    a += 0.5 * ( y1 + y2 ) * ( x1 - x2 ) ;
                    //for (int j(0) ; j<Bodies[Monitored[i][2]].nb_nodes ; j++)
                    //{
                    //    a += Bodies[Monitored[i][2]].nodes[j].x_mass * Bodies[Monitored[i][2]].nodes[j].jacobian ;
                    //}
                    //a = a / Bodies[Monitored[i][2]].density ;
                    current_monitoring.push_back(a) ;
                }
                else if (Monitored[i][1] == 7)
                {
                    current_monitoring.push_back(Bodies[Monitored[i][2]].temperature) ;
                }
            }
            else if (Monitored[i][0] == 15)
            {
                //cout << Monitored[i][0] << Monitored[i][1] << Monitored[i][2] << Monitored[i][3] << Monitored[i][4] << endl ;
                if (Monitored[i][1] == 0)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].gapn ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 1)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].gapt ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 2)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].xsi ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 3)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].xnorm ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 4)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].ynorm ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 5)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].internal[0] ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 6)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].fx ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 7)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].fy ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 8)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].length ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 9)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].contact_elements[j].bodyM ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 10)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = k ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 11)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].nodes[Bodies[k].contact_elements[j].nodeS].x_current ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
                else if (Monitored[i][1] == 12)
                {
                    double q = 0. ;
                    for (int k(0) ; k<Nb_bodies ; k++)
                    {
                        for (int j(0) ; j<Bodies[k].nb_contact_elements ; j++)
                        {
                            if ((Bodies[k].contact_elements[j].fx == 0.) && (Bodies[k].contact_elements[j].fy == 0.))
                                continue ;
                            q = Bodies[k].nodes[Bodies[k].contact_elements[j].nodeS].y_current ;
                            current_monitoring.push_back(q) ;
                        }
                    }
                }
            }
        }
        spying[n].push_back(current_monitoring) ;
    }
}

#endif
