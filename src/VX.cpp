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

    int    loadcase; 
    double tcycle;
    double tequi;
    double totalTime;
    double deltat_max;
    double pull_begin;
    double pull_stretchduration;
    double push_stretchduration;
    double push_begin;
    double push_mag;
    double pull_mag;
    bool applybasalTractionFlag;
    bool perturbseedbucklingflag;
    bool periodicboxflag;

    // Applied deformation
    double DelF00;
    double DelF01;
    double DelF10;
    double DelF11;

    // Re-stiffening parameters
    double memb_thresh_high;
    double memb_thresh_low;
    double memb_lam_high;
    double memb_lam_low;
    int    MAXITER;
    int    nIterChange;
    double SOLTOL;
    double RESTOL;

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
        config.readInto(periodicboxflag,         "periodicboxflag");
        config.readInto(perturbseedbucklingflag, "perturbseedbucklingflag");

        // Applied deformation
        config.readInto(DelF00, "DelF00");
        config.readInto(DelF01, "DelF01");
        config.readInto(DelF10, "DelF10");
        config.readInto(DelF11, "DelF11");

        // Re-stiffening parameters
        config.readInto(memb_thresh_high, "memb_thresh_high");
        config.readInto(memb_thresh_low, "memb_thresh_low");
        config.readInto(memb_lam_high, "memb_lam_high");
        config.readInto(memb_lam_low, "memb_lam_low");

        // Simulation parameters
        config.readInto(loadcase, "loadcase");
        config.readInto(tequi, "tequi");
        config.readInto(tcycle, "tcycle");
        config.readInto(totalTime, "totalTime");
        config.readInto(deltat_max, "deltat_max");
        config.readInto(pull_begin, "pull_begin");
        config.readInto(pull_stretchduration, "pull_stretchduration");
        config.readInto(push_stretchduration, "push_stretchduration");
        config.readInto(push_begin, "push_begin");
        config.readInto(push_mag, "push_mag");
        config.readInto(pull_mag, "pull_mag");

        // Solver parameters
        config.readInto(MAXITER, "MAXITER");
        config.readInto(SOLTOL, "SOLTOL");
        config.readInto(RESTOL, "RESTOL");
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
    double pull_end        = pull_begin + pull_stretchduration;
    double push_end        = push_begin + push_stretchduration;

    // Imposed total deformation gradient on the boundary
    vector<double> Id     = {1.0, 0.0, 0.0, 1.0}; // Identity
    vector<double> Fapp   = {0.0, 0.0, 0.0, 0.0}; // Stores the total imposed deformation on the tissue
    vector<double> Fincr  = {DelF00, DelF01, DelF10, DelF11}; // total increment in F to be imposed

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
        double incr_begin = Getloadstep(tSimuc,pull_begin,pull_end,push_begin,push_end,loadcase,push_mag,pull_mag);
        double incr_end   = Getloadstep(tSimucpdeltatc,pull_begin,pull_end,push_begin,push_end,loadcase,push_mag,pull_mag);
        loadIncr = (incr_end - incr_begin);

        // Update previous timestep solution
        for (int i = 0; i < tissuemesh.nItem(); i++)
        {
            tissuemesh.fields[2 * i + 0]->nodeDOFs0->setValue(tissuemesh.fields[2 * i + 0]->nodeDOFs);
            tissuemesh.fields[2 * i + 1]->nodeDOFs0->setValue(tissuemesh.fields[2 * i + 1]->nodeDOFs);
        }
        
        // ---------------------------------------------------------------------------
        // ----------- Incrementing boundary conditions on dirichlet nodes -----------
        // ---------------------------------------------------------------------------

	if(tissuemesh.myRank() == 0)
		cout << "Applying dirichlet BC..." << endl;

        // Store the master increments in the incrDOFs
        auto incrDOFs = hiperProbl->linProbl->GetLHS();
	incrDOFs->PutScalar(0.0);
        int  hpIdx    = 0;

        for (int dofshandID = 0; dofshandID < hiperProbl->_dhands.size(); dofshandID++) // going through all dofHandlers (DOF and gDOF)
        {
            auto dofshand = hiperProbl->_dhands[dofshandID];

            for (int j = 0; j < dofshand->mesh->loc_nPts(); j++) // getting the number of points in that mechanics (which is equivalent the number of point in the mesh)
            {
                for (int dof = 0; dof < dofshand->numDOFs(); dof++) // getting to the degree of freedom of that dofHandlerS (3 for DOF and more for gDOF)
                {
                    if (dofshandID % 2 == 0)
                    {
                        //~ cout << "mecID : " << mecID << ", j : " << j << ", dof " << dof << endl;
                        int meshID = dofshandID/2;
                        int meshIDloc = tissuemesh.locIdx(meshID);
                        int loc_nPts  = tissuemesh.meshes[meshIDloc].nPts;

                        if (tissuemesh.isItemInPart(meshID))
                        {
                            // Get X
                            vector<double> X = {0.0, 0.0, 0.0};
                            X[0] = tissuemesh.fields[dofshandID]->mesh->nodeCoord(j, 0, IndexType::Global);
                            X[1] = tissuemesh.fields[dofshandID]->mesh->nodeCoord(j, 1, IndexType::Global);
                            X[2] = tissuemesh.fields[dofshandID]->mesh->nodeCoord(j, 2, IndexType::Global);

                            // get increment value
                            vector<double> x_incr = {0.0, 0.0, 0.0};
                            x_incr[0] = loadIncr * (Fincr[0] * X[0] + Fincr[1] * X[1]);
                            x_incr[1] = loadIncr * (Fincr[2] * X[0] + Fincr[3] * X[1]);
                            x_incr[2] = 0.0;

                            // Apply the affine deformation if the node is flagged as dirichlet
                            if (tissuemesh.meshes[meshIDloc].cn_nodes[loc_nPts * dof + j] > 0)
                            {
                                (*incrDOFs)[0][hpIdx] = x_incr[dof];
                            }
                            else
                            {
                                (*incrDOFs)[0][hpIdx] = 0.0;
                            }
                        }
                    }
                    else
                    {
                        int meshID = dofshandID/2;
                        if (tissuemesh.isItemInPart(meshID))
                            (*incrDOFs)[0][hpIdx] = 0.0;
                    }
                    hpIdx++;
                }
            }
        }

        // Update nodeDOFs of slaves nodes following the linear constraints involving dirichlet masters
        hiperProbl->linProbl->SetLHS(incrDOFs);
        hiperProbl->AddSolutionToDOFs(1.0);
        
        // ---------------------------------------------------------------------------
        // --------- Incrementing boundary conditions on periodic simulation ---------
        // ---------------------------------------------------------------------------

        // Apply affine deformation on all nodes resulting in increase in the box size
        if (periodicboxflag == 1)
        {
            if(tissuemesh.myRank() == 0)
                cout << "Applying periodic BC..." << endl;
            
                if (tissuemesh.myRank() == 0)
                cout << "**** NOTE: pulling periodic box with load increment: " << loadIncr << endl;

            for (int i = 0; i < tissuemesh.nItem(); i++)
            {
                // Increment boundary conditions
                if (tissuemesh.isItemInPart(i))
                {
                    int locI = tissuemesh.locIdx(i);
                    int loc_nPts = tissuemesh.meshes[locI].nPts;

                    for (int j = 0; j < loc_nPts; j++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            vector<double> x = {0.0, 0.0, 0.0};
                            x[0] = tissuemesh.fields[2 * i]->nodeDOFs->getValue(0, j, IndexType::Global);
                            x[1] = tissuemesh.fields[2 * i]->nodeDOFs->getValue(1, j, IndexType::Global);
                            x[2] = tissuemesh.fields[2 * i]->nodeDOFs->getValue(2, j, IndexType::Global);

                            vector<double> X = {0.0, 0.0, 0.0};
                            X[0] = tissuemesh.fields[2 * i]->mesh->nodeCoord(j, 0, IndexType::Global);
                            X[1] = tissuemesh.fields[2 * i]->mesh->nodeCoord(j, 1, IndexType::Global);
                            X[2] = tissuemesh.fields[2 * i]->mesh->nodeCoord(j, 2, IndexType::Global);

                            // Update increments
                            vector<double> x_incr = {0.0, 0.0, 0.0};
                            x_incr[0] = loadIncr * (Fincr[0] * X[0] + Fincr[1] * X[1]);
                            x_incr[1] = loadIncr * (Fincr[2] * X[0] + Fincr[3] * X[1]);
                            x_incr[2] = 0.0;

                            tissuemesh.fields[2 * i]->nodeDOFs->setValue(n, j, IndexType::Local, x[n] + x_incr[n]);
                        }
                    }
                }
            }
        }
// Perturb the initial seed for unconstrained nodes to access first buckling mode
int coutSeed = 0;
if (perturbseedbucklingflag == 1)
        {
            if (simTime < 1.0)
                paramStr->dparam[9] = 0.0001;
            else
                paramStr->dparam[9] = eta_f;


            for (int i = 0; i < tissuemesh.nItem(); i++)
            {
                if (tissuemesh.isItemInPart(i))
                {
                    int locI = tissuemesh.locIdx(i);
                    int loc_nPts = tissuemesh.meshes[locI].nPts;

                    for (int j = 0; j < loc_nPts; j++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            if (tissuemesh.meshes[locI].cn_nodes[loc_nPts * n + j] == 0)
                            {
                                vector<double> x = {0.0, 0.0, 0.0};
                                x[0] = tissuemesh.fields[2 * i]->nodeDOFs->getValue(0, j, IndexType::Global);
                                x[1] = tissuemesh.fields[2 * i]->nodeDOFs->getValue(1, j, IndexType::Global);
                                x[2] = tissuemesh.fields[2 * i]->nodeDOFs->getValue(2, j, IndexType::Global);

                                vector<double> X = {0.0, 0.0, 0.0};
                                X[0] = tissuemesh.fields[2 * i]->mesh->nodeCoord(j, 0, IndexType::Global);
                                X[1] = tissuemesh.fields[2 * i]->mesh->nodeCoord(j, 1, IndexType::Global);
                                X[2] = tissuemesh.fields[2 * i]->mesh->nodeCoord(j, 2, IndexType::Global);

                                // Set for 75X10 tissue case
                                if (n == 2 && loadIncr < 0.0 && loadVal > 0.65 && abs(X[0]) < 30.0)
                                {
			
                        	    coutSeed = 1;
                                
				    tissuemesh.fields[2 * i]->nodeDOFs->setValue(n, j, IndexType::Local, x[n] + abs(abs(X[0]) - 30.0) / 100.0);
                                }
                            }
                        }
                    }
                }
            }
        }

	if (coutSeed == 1 && tissuemesh.myRank() == 0)
             cout << "PERTURBING SEED" << endl;

        hiperProbl->UpdateGhosts();
        hiperProbl->FillLinearSystem();
        timerstep.end("Prep");

		// Update simulation status
        if( tissuemesh.myRank() == 0 )
        {
            cout << "Time " << simTime << " of " << totalTime << " with deltat=" << deltat << endl;
            cout << "Starting Newton-Raphson iteration" << endl;
            cout << "Imposed displacement (Fxx): " << 1 + loadVal*Fincr[0] << endl;
            cout << "Friction coefficient : " << paramStr->dparam[9] << endl;
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
