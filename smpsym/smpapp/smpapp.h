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

#ifndef SMPAPP_H
#define SMPAPP_H

#include "smp_addons.h"
#include "smp_simulator.h"
#include "specimen.h"
#include "bracket_model.h"
#include "smp_imgui_handler.h"
#include "smpup_visual.h"
#include "smp_simulation_handler.h"

#ifdef USE_GEOGRAM
#include "smp_remesh_handler.h"
#endif

#ifdef USE_MATLAB
#include "utils/smp_matlab.h"
#endif

#include "cinder/Camera.h"
#include "cinder/GeomIo.h"
#include "cinder/ImageIo.h"
#include "cinder/CameraUi.h"
#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/params/Params.h"
#include "cinder/Log.h"
#include "cinder/app/Platform.h"

#include "cinder/CinderImGui.h"

#include "cinder/Xml.h"

using namespace ci;
using namespace ci::app;
using namespace std;
using namespace smpup;

// Smp application
class SmpApp : public App, SmpAppVisual, private SmpGuiHandler
{
public:
	SmpApp()
#ifdef USE_GEOGRAM
		: mSmpRemeshHandler(this)
#endif
	{ }

	void setup() override;
	void resize() override;
	void update() override;
	void draw() override;

	void mouseDown(MouseEvent event) override;
	void mouseUp(MouseEvent event) override;
	void mouseDrag(MouseEvent event) override;
	void mouseMove(MouseEvent event) override;
	void mouseWheel(MouseEvent event) override;
	void keyDown(KeyEvent event) override;

	void fileDrop(FileDropEvent event) override;

private:

	// Visualization, mouse and keyboard

	Camera*             mCameraActive;
	CameraPersp			mCameraPersp;
	CameraOrtho         mCameraOrtho;
	ivec2               mCameraOrthoMousePos;
	double              mOrthoWidth;
	CameraUi			mCamUi;
	bool				mRecenterCamera;
	vec3				mCameraTarget, mCameraLerpTarget, mCameraViewDirection;
	bool                mRotation;

	gl::BatchRef		mActuatorsBatch;
	int                 mPickedActuator;

	ivec2				mMousePos;
	bool                mMouseLeftButtonDown;  // mouse left button is pressed
	bool                mMouseRightButtonDown; // mouse right button is pressed
	bool                mMouseLeftClickDown;   // mouse left button is pressed, no drag
	bool                mMouseRightClickDown;  // mouse right button is pressed, no drag
	double				mLastMouseLeftDownTime;
	double				mLastMouseRightDownTime;

	void update_ActuatorsBatch();

	// Smpup

	SmpOptions                           mSmpOptions;
	std::shared_ptr<LinearMemModel>      mLinearMemModel;
	std::shared_ptr<BracketModel>        mBracket;
	std::shared_ptr<BracketConfigurer>   mBracketConfigurer;
	std::shared_ptr<CoreSettings>        mCoreSettings;
	float                                mTimestep;

#ifdef USE_GEOGRAM
	SmpRemeshHandler mSmpRemeshHandler;
#endif

	std::list<std::shared_ptr<SimulationHandler>>   mSimulationHandlers;
	std::shared_ptr<SimulationHandler>              mSelectedSimHandler;

	// Remove all simulations and clear consol
	void remove_all();

	// Import surface to generate stencils
	void import_surface(fs::path path = fs::path());

	// Read stencil mesh info
	void read_file_info(fs::path path = fs::path());

	// Open current stencil mesh and prepare a simulator
	void open_current_file(SmpSettings settings);

	// Create regular stencil mesh and prepare a simulator
	void new_regular(int rows, int cols, double length, double thickness, SmpSettings settings);

	// Create single actuator stencil mesh and prepare a simulator
	void new_single(double length, double thickness, SmpSettings settings);

	// Perofrom picking action for mouse left click
	template<typename DerivedA, typename DerivedB>
	Eigen::Vector3d pick_triangle(const Eigen::DenseBase<DerivedA>& V, const Eigen::DenseBase<DerivedB>& F, int& picked);

	// Perofrom picking action for mouse left click
	int performActuatorPicking();

	int performBasePicking();

	// Add simulation handler
	int addSimulationHandler(std::string name, SmpSettings settings, const Smp& smp);

	// Remove simulation handler
	int removeSimulationHandler(std::shared_ptr<SimulationHandler>& shandler);


	// ImGui and controls

	bool                mPlayback;
	double              mPlaySpeed;
	bool                mPlayLoop;
	bool                mLockRestConfiguration;
	bool                mPlaySmooth;

	int                 mCurrentFrame;
	double              mCurrentFrameAlpha; // interpolate between frames for smoothness

	bool                mExitRequest;

	std::string         mCurrentFile;

	int                 mTimeLandscapeMax;
	int                 mTimeLandscapeMin;
	enum                TimeBrushMode { SCULPT, SMOOTHEN };
	TimeBrushMode       mTimeBrushMode;
	int                 mTimeLandscapeRadius;
	int                 mTimeLandscapeIncrement;

	std::stringstream   buf_console;

#ifdef USE_MATLAB
public:
	bool                matlab_started; // is true once Matlab engine has been started
	Matlab              matlab;
#endif

private:

	//Gui popups
	void PopupClearAll();
	void PopupNewSingle();
	void PopupNewRegular();
	void PopupOpenStencil();
	void PopupCloseApp();
	void PopupBuildMesh();
	void PopupBuildVideo();
	void PopupCoreSettings();
	void PopupMembraneLinearity();
	void PopupBuildSpecimenBlock();
	void PopupBracketSampling();

	// Gui windows
	void GuiMainWindow();
	void GuiSurfaceImport();
	void GuiLandscapeEditor();
	void GuiSolverInfo();
	void GuiConsole();
	void GuiBracketExplorer();
	void GuiPlottingWindow();
	void GuiDataExportWindow();
	void GuiDebugWindow();

	void GuiMainMenu();
	void GuiPopups();

	void showGui();

#if ! defined( CINDER_GL_ES )
	params::InterfaceGlRef	mParams;

	typedef std::vector<params::InterfaceGl::OptionsBase> ParamGroup;
	std::vector<ParamGroup>	mPrimitiveParams;
#endif
};

#endif
