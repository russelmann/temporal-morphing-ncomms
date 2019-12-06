// This file is part of smpsym, an inverse design and simulation tool for
// research paper Guseinov R. et al "Programming temporal morphing of
// self-actuated shells"
//
// Copyright (C) 2019 Ruslan Guseinov <guseynov.ruslan@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "smpapp.h"
#include "smp_imgui.h"

void SmpApp::GuiMainMenu()
{
	if (!ui::BeginMainMenuBar()) return;

	if (ui::BeginMenu("File"))
	{
		if (ui::MenuItem("New single...")) modal_new_single_request = true;
		if (ui::MenuItem("New regular...")) modal_new_regular_request = true;
		if (ui::MenuItem("Open...", "Ctrl+O")) { read_file_info(); model_open_stencil_request = true; }
		ui::Spacing();
		if (mSelectedSimHandler && ui::MenuItem("Save simulation..."))
		{
			fs::path path = app::getSaveFilePath(mSelectedSimHandler->simulator->get_smp()->working_dir + "/" + mSelectedSimHandler->simulator->get_name() + ".bin",
				vector<string>{"json", "bin"});
			if (!path.empty()) mSelectedSimHandler->simulator->save(path.string());
		}
		if (ui::MenuItem("Load simulation..."))
		{
			fs::path path = app::getOpenFilePath(fs::path(), vector<string>{"json", "bin"});
			if (!path.empty())
			{
				auto simulator = make_unique<SmpSimulator>(mLinearMemModel, mBracket);
				if (!simulator->load(path.string()))
				{
					clog << "Failed to read file\n";
					return;
				}
				simulator->set_smp_core_settings(mCoreSettings);
				simulator->set_smp_working_dir(path.parent_path().string());

				shared_ptr<SimulationHandler> shandler = make_shared<SimulationHandler>((SmpAppVisual*)this);

				shandler->visible = true;

				shandler->color = Color(simulator->get_settings().color(0), simulator->get_settings().color(1), simulator->get_settings().color(2));

				shandler->simulator.reset(simulator.release());
				shandler->time_landscape = shandler->simulator->get_smp()->time_landscape;

				shandler->createStencils();
			    shandler->update_MembraneBatch(mCurrentFrame, mCurrentFrameAlpha);
				shandler->update_RigidBodyBatch(mCurrentFrame, mCurrentFrameAlpha);
				shandler->update_TimeLandscapeBatch();
				update_ActuatorsBatch();

				mSimulationHandlers.emplace_front(shandler);
				mSelectedSimHandler = mSimulationHandlers.front();

				mPickedActuator = -1;
				clog << "File ... is successully loaded\n";
			}
		}
		ui::Spacing();
		if (ui::MenuItem("Remove all simulations")) modal_clear_all_request = true;
		if (ui::MenuItem("Exit", "Alt+X")) mExitRequest = true;
		ui::EndMenu();
	}

	if (mSelectedSimHandler && ui::BeginMenu("Export"))
	{
		if (ui::MenuItem("Build mesh...")) modal_build_mesh_request = true;
		if (ui::MenuItem("Build video...")) modal_build_video_request = true;
		ui::EndMenu();
	}

	if (ui::BeginMenu("View"))
	{
		ui::Checkbox("Surface import", &mShowSurfaceImport);
		ui::Spacing();
		ui::Checkbox("Time landscape editor", &mShowTimeLandscapeEditor);
		ui::Checkbox("Solver info", &mShowCurrentSolverInfo);
		ui::Spacing();
		ui::Checkbox("Bracket explorer", &mShowBracketExplorer);
		ui::Checkbox("Plotting", &mShowPlottingWindow);
		ui::Checkbox("Data export", &mShowDataExportWindow);
		ui::Checkbox("Console", &mShowConsole);
		ui::Spacing();
		ui::Checkbox("Debug window", &mShowDebugWindow);
		ui::EndMenu();
	}

	if (ui::BeginMenu("Tools"))
	{
		if (ui::MenuItem("Membrane linearity sampling...")) modal_membrane_linearity_request = true;
		if (ui::MenuItem("Build specimen block...")) modal_build_specimen_block_request = true;
		if (ui::MenuItem("Bracket sampling...")) modal_bracket_sampling_request = true;

		if (ui::MenuItem("Core settings...")) modal_core_settings_request = true;

		ui::EndMenu();
	}

#ifdef USE_MATLAB
	if (ui::BeginMenu("Matlab"))
	{
		if (!matlab_started) ui::Text("Waiting for matlab engine...");
		else
		{
			if (!mSelectedSimHandler) ui::Text("No simulation");
			else
			{
				if (ui::MenuItem("Plot actuator length histogram"))
				{
					matlab.figure("act_len", "Actuator initial lengths", "length, mm", "count", false);
					matlab.histogram(mSelectedSimHandler->simulator->get_actuator_lengths_at_frame(mCurrentFrame), 20, Eigen::Vector3d{ 0, 0, 1 });
				}

				if (mPickedActuator >= 0)
				{
					ui::Separator();
					ui::Text("For selected actuator:");

					if (ui::MenuItem("   Plot length vs time"))
					{
						Eigen::VectorXd times = (natural_seriesd(mSelectedSimHandler->simulator->get_simulated_frames()).array() - 1) * mSelectedSimHandler->simulator->get_timestep();
						matlab.figure("actuator_length", "Actuator length", "time, sec", "length, mm", false);
						matlab.plot(times, mSelectedSimHandler->simulator->get_actuator_lengths(mPickedActuator).head(times.size()));
					}

					if (ui::MenuItem("   Plot rest length vs time"))
					{
						Eigen::VectorXd times = (natural_seriesd(mSelectedSimHandler->simulator->get_simulated_frames()).array() - 1) * mSelectedSimHandler->simulator->get_timestep();
						Eigen::MatrixXd restlen(times.size(), 4);
						static Eigen::Vector4d len;
						for (int i = 0; i < restlen.rows(); ++i)
						{
							mSelectedSimHandler->simulator->get_bracket_restlen_at_frame(mPickedActuator, i, len);
							restlen.row(i) = len.transpose();
						}

						matlab.figure("act_restlen", "Actuator rest lengths", "time, sec", "length, mm");

						for (int i = 0; i < 4; ++i)
							matlab.plot(times, restlen.col(i), i < 2 ? "blue" : "green");
					}

					if (ui::MenuItem("   Plot forces vs time"))
					{
						Eigen::VectorXd times = (natural_seriesd(mSelectedSimHandler->simulator->get_simulated_frames()).array() - 1) * mSelectedSimHandler->simulator->get_timestep();
						Eigen::MatrixXd forces(times.size(), 4);
						static Eigen::Vector4d bracket_forces;
						for (int i = 0; i < forces.rows(); ++i)
						{
							mSelectedSimHandler->simulator->get_bracket_forces_at_frame(mPickedActuator, i, bracket_forces);
							forces.row(i) = bracket_forces.transpose();
						}

						matlab.figure("act_force", "Actuator forces", "time, sec", "force, N");

						for (int i = 0; i < 4; ++i)
							matlab.plot(times, forces.col(i), i < 2 ? "blue" : "green");
					}

					if (ui::MenuItem("   Plot stress-strain"))
					{
						Eigen::VectorXd lengths(mSelectedSimHandler->simulator->get_simulated_frames());
						Eigen::VectorXd forces(lengths.size());
						for (int i = 0; i < forces.rows(); ++i)
						{
							Eigen::Vector4d bracket_forces;
							mSelectedSimHandler->simulator->get_bracket_forces_at_frame(mPickedActuator, i, bracket_forces);
							forces(i) = bracket_forces.sum();
							lengths(i) = mSelectedSimHandler->simulator->get_actuator_length_at_frame(mPickedActuator, i);
						}

						matlab.figure("act_stress_strain", "Actuator membrane stress-strain", "length, mm", "force, N");

						if (mCurrentFrame > 0)
						{
							// remove first values since they represent currently rest state
							lengths = Eigen::VectorXd(lengths.bottomRows(lengths.rows() - 1));
							forces = Eigen::VectorXd(forces.bottomRows(forces.rows() - 1));

							// remove curve after actuation
							// TODO: use simulator->collided
							lengths = Eigen::VectorXd(lengths.topRows(mCurrentFrame - 1));
							forces = Eigen::VectorXd(forces.topRows(mCurrentFrame - 1));
						}

						matlab.plot(lengths, forces, "blue");

						// Output matrix for stress-strain membrane mimicing springs fitting
						//Eigen::MatrixXd outm(lengths.rows(), 3);
						//outm.col(0) = lengths;
						//outm.col(1).setConstant((mSelectedSimHandler->simulator->finfo()[0]->body_origins.row(mSelectedSimHandler->smp->FF(mPickedActuator, 0)) -
						//	                     mSelectedSimHandler->simulator->finfo()[0]->body_origins.row(mSelectedSimHandler->smp->FF(mPickedActuator, 1))).norm());
						//outm.col(2) = forces;
					} tooltip("Stress-strain curve for total force acting on brackets vs current displacement");
				}
			}
		}

		ui::EndMenu();
	}
#endif
	ui::EndMainMenuBar();
}
