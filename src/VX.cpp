/*
    *******************************************************************************
    Copyright (c) 2017-2025
    Authors/Contributors: Adam Ouzeri, Sohan Kale, Alejandro Torres-Sánchez, Daniel Santos-Oliván
    *******************************************************************************
    This file is part of agvm
    Project homepage: https://github.com/aouzeri/multiscale-epithelial-mechanics-PRXLife
    Distributed under the GNU General Public License, see the accompanying
    file LICENSE or https://opensource.org/license/gpl-3-0.
    *******************************************************************************
*/

#include <iostream>
#include <mpi.h>

#include "hl_DistributedClass.h"
#include "hl_TypeDefs.h"
#include "hl_DistributedMesh.h"
#include "hl_ParamStructure.h"
#include "hl_LoadMesh.h"
#include "hl_Math.h"
#include "hl_ConfigFile.h"
#include "hl_LocMongeParam.h"
#include "hl_StructMeshGenerator.h"
#include "hl_MeshLoader.h"

#include "hl_ConsistencyCheck.h"
#include <dmumps_c.h>
#include "hl_LinearSolver_Direct_MUMPS.h"
#include "hl_NonlinearSolver_NewtonRaphson.h"

#include "hl_DOFsHandler.h"
#include "AuxVX.h"

#include <random>
#include <hl_Timer.h>
#include "hl_HiPerProblem.h"

using namespace std;
using namespace hiperlife;

int main(int argc, char *argv[])
{
    hiperlife::Init(argc, argv);
    const int myRank = hiperlife::MyRank();

    // ---------------------------------------------------------------------------
    // ----------------------- Setting simulation parameters ---------------------
    // ---------------------------------------------------------------------------

    // Initializing set of parameters
    int    gPts;
    double deltat;
    int nCells;
    double kp;
    double kd;
    double ConcInit;
    double bvisc;
    double lame1;
    double lame2;
    double activA;
    double activB;
    double activL;
    double eta_f;
    double basaltractions;
    double tcycle;
    double tequi;
    double inflation_begin;
    double inflation_duration;
    double deflation_duration;
    double deflation_begin;
    double totalTime;
    double deltat_max;
    double deflation_mag;
    double inflation_mag;
    bool applybasalTractionFlag;

    // Re-stiffening parameters
    double memb_thresh_high;
    double memb_thresh_low;
    double memb_lam_high;
    double memb_lam_low;

    int    MAXITER;
    int    nIterChange;

    double SOLTOL;
    double RESTOL;
    int    nResSteps;

    // Model parameters and solver options are read from a configuration file
    if (argc > 1)
    {
        if( myRank == 0 )
            cout << "Reading the config file: " << argv[1] << endl;
        const char *config_filename = argv[1];
        ConfigFile config(config_filename);

        config.readInto(gPts, "gPts");
        config.readInto(deltat, "deltat");
        config.readInto(nCells, "nCells");

        // Material parameters
        config.readInto(kp, "kp");
        config.readInto(kd, "kd");
        config.readInto(ConcInit, "ConcInit");
        config.readInto(bvisc, "bvisc");
        config.readInto(lame1, "lame1");
        config.readInto(lame2, "lame2");
        config.readInto(activA, "activA");
        config.readInto(activB, "activB");
        config.readInto(activL, "activL");
        config.readInto(eta_f, "eta_f");
        config.readInto(basaltractions, "basaltractions");
        config.readInto(applybasalTractionFlag, "applybasalTractionFlag");

        // Re-stiffening parameters
        config.readInto(memb_thresh_high, "memb_thresh_high");
        config.readInto(memb_thresh_low, "memb_thresh_low");
        config.readInto(memb_lam_high, "memb_lam_high");
        config.readInto(memb_lam_low, "memb_lam_low");

        // Simulation parameters
        config.readInto(tequi, "tequi");
        config.readInto(tcycle, "tcycle");
        config.readInto(inflation_begin, "inflation_begin");
        config.readInto(inflation_duration, "inflation_duration");
        config.readInto(deflation_duration, "deflation_duration");
        config.readInto(deflation_begin, "deflation_begin");
        config.readInto(totalTime, "totalTime");
        config.readInto(deltat_max, "deltat_max");
        config.readInto(deflation_mag, "deflation_mag");
        config.readInto(inflation_mag, "inflation_mag");

        // Solver parameters
        config.readInto(MAXITER, "MAXITER");
        config.readInto(SOLTOL, "SOLTOL");
        config.readInto(RESTOL, "RESTOL");
        config.readInto(nResSteps, "nResSteps");
        config.readInto(nIterChange, "nIterChange");

    }
    else
    {
        if( myRank == 0 )
            cout << "No config file provided, please provide a config file..." << endl;
        throw runtime_error("Please provide a config file in the run command");
    }

    // ---------------------------------------------------------------------------
    // ----------- Read mesh from input files and generate tissuemesh ------------
    // ---------------------------------------------------------------------------

    Timer timer({"total"});
    timer.start("total");

    if( myRank == 0 )
        cout << "Generating tissuemesh..." << endl;
    tissueMesh tissuemesh;


    if( myRank == 0)
	cout << "Loading mesh..." << endl;
    tissuemesh.loadMesh("vertexmesh");
  
    if( myRank == 0)
	cout << "Updating tissue mesh..." << endl;
    tissuemesh.Update();
    
    if( myRank == 0 )
        cout << "Finished generating tissuemesh..." << endl;


    // ---------------------------------------------------------------------------
    // --------------- Generate user structure and HiPerProblem ----------------
    // ---------------------------------------------------------------------------

    // Assign parameters in a user structure
    SmartPtr<ParamStructure> paramStr = Create<ParamStructure>();
    paramStr->dparam.resize(16);

    paramStr->dparam[0]  = deltat;
    paramStr->dparam[1]  = kp;
    paramStr->dparam[2]  = kd;
    paramStr->dparam[3]  = bvisc;
    paramStr->dparam[4]  = lame1;
    paramStr->dparam[5]  = lame2;
    paramStr->dparam[6]  = activA; 
    paramStr->dparam[7]  = activB; 
    paramStr->dparam[8] = activL; 
    paramStr->dparam[9] = eta_f;    
    paramStr->dparam[10] = basaltractions; // Updated in time loop 
    paramStr->dparam[11] = ConcInit;
    paramStr->dparam[12] = memb_thresh_high;
    paramStr->dparam[13] = memb_thresh_low;
    paramStr->dparam[14] = memb_lam_high;
    paramStr->dparam[15] = memb_lam_low;


    // Setting simualtion parameters
    int    nSave    = 0;
    double simTime  = 0.0;
    double loadVal  = 0.0;

    double tSimuc, tSimucpdeltatc, loadIncr;

    // Write initial loaded mesh
    if( tissuemesh.myRank() == 0 )
        cout << "writing initial output file..." << endl;
    tissuemesh.saveSolution("solution"+to_string(nSave), paramStr, simTime);      
    nSave ++;  

    if( tissuemesh.myRank() == 0 )
        cout << "Finished writing initial output file..." << endl;

    if( tissuemesh.myRank() == 0 )
        cout << "Creating the hiperproblem..." << endl;	

    SmartPtr<HiPerProblem> hiperProbl   = tissuemesh.generateHiPerProblem(paramStr,gPts);
    
    if( tissuemesh.myRank() == 0 )
        cout << "Finished creating the hiperproblem..." << endl;


   if (tissuemesh.myRank() == 0 )
	  cout << "Applying constraints..."<< endl; 

    // Updating initial conditions
    for (int i = 0; i < tissuemesh.nItem(); i++)
    {
        if (tissuemesh.isItemInPart(i))
        {
            tissuemesh.fields[2 * i + 1]->nodeDOFs->setValue(4, 0, IndexType::Local, ConcInit);
        }
    }

    //Set parameter
    double stepFactor = 0.90;
    bool startConstrainingMesh = false;

    //double parameters for NR
    std::vector<double> dparamNR;
    dparamNR.resize(4);
    dparamNR[0] = SOLTOL;
    dparamNR[1] = RESTOL;

    //Integer parameters for NR
    std::vector<int> iparamNR;
    iparamNR.resize(1);
    iparamNR[0] = MAXITER;

    double inflation_end        = inflation_begin + inflation_duration;
    double deflation_end        = deflation_begin + deflation_duration;

    if (tissuemesh.myRank() == 0 )
          cout << "Initial fill of hiperproblem..."<< endl;

    // Initial Fill
    hiperProbl->UpdateGhosts();
    hiperProbl->FillLinearSystem();

     if (tissuemesh.myRank() == 0 )
          cout << "Finished initial fill of hiperproblem..."<< endl;

    //Solvers
    if( tissuemesh.myRank() == 0 )
        cout << "Creating the linear and nonlinear solver..." << endl;

    // Newton-raphson iterator
    SmartPtr<NewtonRaphsonNonlinearSolver>  nonlinSolver = Create<NewtonRaphsonNonlinearSolver>();

    SmartPtr<MUMPSDirectLinearSolver> linSolver = Create<MUMPSDirectLinearSolver>();
    if( tissuemesh.myRank() == 0 )
        cout << "Using Mumps direct solver..." << endl;
    linSolver->setHiPerProblem(hiperProbl);
    linSolver->setDefaultParameters();
    linSolver->Update();

    nonlinSolver->setLinearSolver(linSolver);
    nonlinSolver->setConvRelTolerance(false);
    nonlinSolver->setMaxNumIterations(MAXITER);
    nonlinSolver->setResTolerance(RESTOL);
    nonlinSolver->setSolTolerance(SOLTOL);
    nonlinSolver->setPrintIntermInfo(true);
    nonlinSolver->setPrintSummary(true);
    nonlinSolver->Update();
    
    if( tissuemesh.myRank() == 0 )
        cout << "Finished creating the solver..." << endl;

    if( tissuemesh.myRank() == 0 )
        cout << "Starting simulation loop..." << endl;

    bool resetDeltat = true;

    ////////////////////////////////////////////////////////////// Time -Loop ////////////////////////////////////////////////////////////// 
    while ( simTime < totalTime )
    {	 

        Timer timerstep({"Step","Prep","Solver","PostProc"});
        timerstep.start("Step");   

        timerstep.start("Prep");
        
        //////////////////////////////////////// Ad hoc modification ////////////////////////////////////////

        // Applying tractions on basal boundary
        if (applybasalTractionFlag == true)
            paramStr->dparam[10] = basaltractions*loadVal;
        else
            paramStr->dparam[10] = 0.0;
    
		//**************************************************************************************************************//

        // Step time inside the cycle
        tSimuc         = fmod(simTime - tequi, tcycle);
        tSimucpdeltatc = fmod(simTime + deltat - tequi, tcycle);

        // Get the load increment to be imposed
        double incr_begin = Getloadstep(tSimuc,inflation_begin,inflation_end,deflation_begin,deflation_end,deflation_mag,inflation_mag);
        double incr_end   = Getloadstep(tSimucpdeltatc,inflation_begin,inflation_end,deflation_begin,deflation_end,deflation_mag,inflation_mag);
        loadIncr = (incr_end - incr_begin);

        // Update previous timestep solution
        for (int i = 0; i < tissuemesh.nItem(); i++)
        {
            tissuemesh.fields[2 * i + 0]->nodeDOFs0->setValue(tissuemesh.fields[2 * i + 0]->nodeDOFs);
            tissuemesh.fields[2 * i + 1]->nodeDOFs0->setValue(tissuemesh.fields[2 * i + 1]->nodeDOFs);
        }     

        hiperProbl->UpdateGhosts();
        hiperProbl->FillLinearSystem();
        timerstep.end("Prep");

		// Update simulation status
        if( tissuemesh.myRank() == 0 )
        {
            cout << "Time " << simTime << " of " << totalTime << " with deltat=" << deltat << endl;
            cout << "Starting Newton-Raphson iteration" << endl;
            cout << "Step number: "<< nSave << endl;
        }

        // Non-linear solver
        timerstep.start("Solver");
        bool ierr = nonlinSolver->solve();
        timerstep.end("Solver");

    
        if (ierr)
        {

		// Update simulation time
            simTime += deltat;

	     // Increment the loadstep value
            loadVal = loadVal + loadIncr;


	    if( tissuemesh.myRank() == 0  )
                cout << "Postprocessing... " << endl;

        timerstep.start("PostProc");

	    if( tissuemesh.myRank() == 0  )
                cout << "Printing outputs... " << endl;


            // Update deltat    
            int nrIter = nonlinSolver->numberOfIterations();
            if ( nrIter < nIterChange)
            {
                deltat /= stepFactor;

                if (deltat > deltat_max)
                {
                    deltat = deltat_max;
                }
            }

            if( tissuemesh.myRank() == 0  )
                cout << "Saving files in time-step " << nSave << " corresponding to solution at time = " << simTime << endl;

            // Save .vtu and .vtm files
            timerstep.start("PostProc");
            tissuemesh.saveSolution("solution"+to_string(nSave), paramStr, simTime);
            nSave ++;  

	   
        }
        else
        {
				
            if ( tissuemesh.myRank() == 0 )
	            std::cout << "ERROR::EXIT from MAIN LOOP with error..." << ierr << endl;

        
            //Modify time-step
	        deltat *= stepFactor;

            // Recover nodeDOFs from nodeDOFs0
            for (int i = 0; i < tissuemesh.nItem(); i++)
            {
                tissuemesh.fields[2*i+0]->nodeDOFs->setValue(tissuemesh.fields[2*i+0]->nodeDOFs0);
                tissuemesh.fields[2*i+1]->nodeDOFs->setValue(tissuemesh.fields[2*i+1]->nodeDOFs0);
                tissuemesh.fields[2*i+0]->nodeDOFs->UpdateGhosts();
                tissuemesh.fields[2*i+1]->nodeDOFs->UpdateGhosts();
            }
            hiperProbl->UpdateGhosts();
        }

        // Update the time-step value
        paramStr->dparam[0] = deltat;

        // Output the timing information
        timerstep.end("Step");
        timerstep.end("PostProc");
        timer.end("total");


        string prefix = "Step time, [Np = " + to_string(tissuemesh.numProcs()) + "] : ";
        timerstep.printTesting(prefix);
        if( tissuemesh.myRank() == 0 )
        {
            cout << "Simulation time (min) : " << timer.printAccumTime("total")/60.0 << endl;
            cout << "-------------"<< endl;
            cout << endl;
        }

       timer.start("total");
    }

    timer.end("total");

    if( tissuemesh.myRank() == 0 )
    {
        cout << "-----------------------------------------------------------------"    << endl;
        cout << "Total simulation time (min) : " << timer.printAccumTime("total")/60.0 << endl;
        cout << endl;
    }

    MPI_Finalize();

    return 0;
}
