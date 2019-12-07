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
#include "CinderImGui.h"

void prepareSettings( App::Settings* settings )
{
	settings->setWindowSize( 1024, 768 );
	settings->setHighDensityDisplayEnabled();
	settings->setMultiTouchEnabled( false );
	settings->disableFrameRate();
}

void SmpApp::remove_all()
{
	buf_console.str(string());
	buf_console.clear();
	mSelectedSimHandler.reset();
	mSimulationHandlers.clear();
}

void SmpApp::setup()
{
	gl::enableVerticalSync( true );

	mMouseLeftButtonDown = false;
	mMouseRightButtonDown = false;
	mMouseLeftClickDown = false;
	mMouseRightClickDown = false;
	mLastMouseLeftDownTime = 0;

	// Initialize variables
	
	filesystem::path asset_dir = app::loadAsset("smp_anchor")->getFilePath().parent_path().c_str();
	mCoreSettings = make_shared<CoreSettings>(asset_dir);
	mCoreSettings->serialization();
	mTimestep = mCoreSettings->timestep;
	mSmpOptions.set_default_directory(mCoreSettings->get_smpup_path().append("model"));

	mExitRequest = false;

	getWindow()->setTitle("SmpApp");

	// Setup the camera.
	mCameraActive = &mCameraPersp;
	mCameraPersp.lookAt(normalize(vec3(3, 3, 6)) * 145.0f, mCameraTarget);
	mCameraOrtho.lookAt(vec3(0, 0, 200), vec3(0, 0, -1));
	mOrthoWidth = 10;
	mCameraOrtho.setOrtho(-mOrthoWidth, mOrthoWidth, -mOrthoWidth / getWindowAspectRatio(), mOrthoWidth / getWindowAspectRatio(), 0.1, 1000);
	mCameraOrtho.setAspectRatio(getWindowAspectRatio());
	
	mCamUi = CameraUi( &mCameraPersp );

	mRotation = false;
	update_grid();

	mPlayback = true;
	mPlaySpeed = 1;
	mPlayLoop = true;
	mLockRestConfiguration = true;
	mPlaySmooth = true;
	mCurrentFrame = 0;
	mCurrentFrameAlpha = 0;
	mPickedActuator = -1;

	mTimeLandscapeMax = 120;
	mTimeLandscapeMin = 30;
	mTimeBrushMode = SCULPT;
	mTimeLandscapeRadius = 50;
	mTimeLandscapeIncrement = 3;

	ui::initialize();

	clog.rdbuf(buf_console.rdbuf());

	// Load linear membrane model
	{
		ifstream ss(mCoreSettings->get_smpup_path().append("model").append("linear_membrane_model.json").c_str());
		if (ss.good())
		{
			cereal::JSONInputArchive archive(ss);
			mLinearMemModel = make_shared<LinearMemModel>();
			archive(cereal::make_nvp("linear_membrane_model", *mLinearMemModel));
		}
		else clog << "Error: Could not load linear membrane model\n";
	}

	// Load bracket model
	{
		ifstream ss(mCoreSettings->get_smpup_path().append("model").append("bracket_model.json").c_str());
		if (ss.good())
		{
			cereal::JSONInputArchive archive(ss);
			mBracket = make_shared<BracketModel>();
			archive(cereal::make_nvp("bracket_model", *mBracket));
		}
		else clog << "Error: Could not load bracket model\n";
	}

	// Load bracket configurer
	mBracketConfigurer = make_shared<BracketConfigurer>(mLinearMemModel, mBracket);
	mBracketConfigurer->set_default_directory(mCoreSettings->get_smpup_path().append("model"));
	if (!mBracketConfigurer->load())
        clog << "Error: Could not load bracket configurer\n";

#ifdef USE_MATLAB
	new thread([this]() {
		matlab_started = false;
        matlab_started = matlab.start();
        if (matlab_started) clog << "Matlab engine started\n";
        else clog << "Error loading Matlab engine\n";
	});
#endif
}

int SmpApp::addSimulationHandler(string name, SmpSettings settings, const Smp& smp)
{
	string cname = name;
	int iname = 1;
	while (true)
	{
		bool exists = false;
		for (auto shandler : mSimulationHandlers)
			if (cname == shandler->simulator->get_name())
			{
				exists = true;
				break;
			}
		if (exists)
        {
            cname = name + ' ' + string(3 - to_string(iname).length(), '0') + to_string(iname);
            ++iname;
        }
		else break;
	}

	if ((settings.color.array() < 0).any() || (settings.color.array() > 1).any())
		settings.color = Eigen::Vector3f::Random().array().abs();

	shared_ptr<SimulationHandler> shandler = make_shared<SimulationHandler>((SmpAppVisual*)this);

	shandler->visible = true;

    shandler->time_landscape.setConstant(smp.V.rows(), 30);
	
	shandler->simulator = make_unique<SmpSimulator>(0, 120, mTimestep, settings, mLinearMemModel, mBracket);
	shandler->simulator->set_name(cname);
	shandler->simulator->load_smp(smp);
	//shandler->simulator.reset(shandler->simulator->configure_brackets(*mBracketConfigurer, shandler->time_landscape).release());
	shandler->color = Color(settings.color(0), settings.color(1), settings.color(2));
	shandler->createStencils();

	mSimulationHandlers.emplace_front(shandler);
	mSelectedSimHandler = shandler;
	clog << "Simulation '" + shandler->simulator->get_name() + "' successfully added\n";

	shandler->update_MembraneBatch(mCurrentFrame, mCurrentFrameAlpha);
	shandler->update_RigidBodyBatch(mCurrentFrame, mCurrentFrameAlpha);
	shandler->update_TimeLandscapeBatch();
	update_ActuatorsBatch();

	return 0;
}

int SmpApp::removeSimulationHandler(shared_ptr<SimulationHandler>& shandler)
{
	if (mSelectedSimHandler == shandler) mSelectedSimHandler.reset();
	string name = shandler->simulator->get_name();
	mSimulationHandlers.remove(shandler);
	shandler.reset();
	clog << "Simulation '" << name << "'  removed\n";
	if (!mSimulationHandlers.empty()) mSelectedSimHandler = mSimulationHandlers.front();
	return 0;
}

void SmpApp::import_surface(fs::path path)
{
	if (path.empty())
	{
		path = app::getOpenFilePath(fs::path(), vector<string>{"obj"});
		if (path.empty()) return;
	}
#ifdef USE_GEOGRAM
	string abs_path = path.string();
	mSmpRemeshHandler.geo_remesher->load_mesh(abs_path);
	mSmpRemeshHandler.update_ImportSurfacesBatch();
#else
	clog << "Import is not possible. This version is compiled without Geogram.\n";
#endif
}

void SmpApp::read_file_info(fs::path path)
{
	if (path.empty())
	{
		path = app::getOpenFilePath(fs::path(), vector<string>{"obj"});
		if (path.empty()) return;
	}
	clog << "File <" << path.filename() << "> is opened\n";
	mCurrentFile = path.string();
}

void SmpApp::open_current_file(SmpSettings settings)
{
	clog << "Opening stencil structure\n";

	mSmpOptions.serialization();

	Smp smp;
	smp.core_settings = const_pointer_cast<const CoreSettings>(mCoreSettings);
	smp.options = mSmpOptions;
	if (!smp.read_mesh(mCurrentFile))
	{
		clog << "Failed to read mesh\n";
		return;
	}
	
	if (addSimulationHandler(smp.name, settings, smp))
	{
		clog << "Error: file cannot be opened\n";
		return;
	}

	mPickedActuator = -1;
	clog << "File <" << mCurrentFile << "> is successully loaded\n";
}

void SmpApp::new_regular(int rows, int cols, double length, double thickness, SmpSettings settings)
{
	clog << "Initializing regular structure\n";

	mSmpOptions.serialization();

	mSmpOptions.min_gap = length / 3; // assuming 1.5 contraction
	mSmpOptions.min_bumper = length / 3 - mSmpOptions.base_r;

	Smp smp;
	smp.core_settings = const_pointer_cast<const CoreSettings>(mCoreSettings);
	smp.options = mSmpOptions;
	smp.init_regular(rows, cols);
	smp.spth.setConstant(thickness);

	if (addSimulationHandler(smp.name, settings, smp))
	{
		clog << "Error: file cannot be opened\n";
		return;
	}

	mPickedActuator = -1;
	clog << "Regular structure is successully loaded\n";
}

void SmpApp::new_single(double length, double thickness, SmpSettings settings)
{
	clog << "Initializing single actuator\n";
	
	mSmpOptions.serialization();

	mSmpOptions.min_gap = length / 3; // assuming 1.5 contraction
	mSmpOptions.min_bumper = length / 3 - mSmpOptions.base_r;

	Smp smp;
	smp.core_settings = const_pointer_cast<const CoreSettings>(mCoreSettings);
	smp.options = mSmpOptions;
	smp.init_single();
	smp.spth.setConstant(thickness);

	if (addSimulationHandler(smp.name, settings, smp))
	{
		clog << "Error: linear membrane initialization failure\n";
		return;
	}

	mPickedActuator = -1;
	clog << "Regular structure is successully loaded\n";
}

void SmpApp::update()
{
	// Controls

	if (mMouseRightButtonDown && mSelectedSimHandler && mShowTimeLandscape)
	{
		static double u_last_rightbutton_update = 0;
		if (u_last_rightbutton_update < getElapsedSeconds() - 0.2)
		{
			const Smp* smp = mSelectedSimHandler->simulator->get_smp();
			Eigen::MatrixXd V;
			V.resize(smp->uv.rows(), 3);
			V.leftCols(2) = smp->uv;
			V.col(2) = mSelectedSimHandler->time_landscape;
			int picked;
			Eigen::Vector3d picked_point = pick_triangle(V, smp->F, picked);
			if (0 <= picked)
			{
				if (mTimeBrushMode == SCULPT)
				{
					for (int i = 0; i < smp->uv.rows(); ++i)
					{
						double dist = (picked_point.head(2).transpose() - smp->uv.row(i)).norm();
						if (mTimeLandscapeRadius < dist) continue;
						mSelectedSimHandler->time_landscape(i) += double(mTimeLandscapeIncrement) * (mTimeLandscapeRadius - dist) / mTimeLandscapeRadius;
						mSelectedSimHandler->time_landscape(i) = min(mSelectedSimHandler->time_landscape(i), double(mTimeLandscapeMax));
						mSelectedSimHandler->time_landscape(i) = max(mSelectedSimHandler->time_landscape(i), double(mTimeLandscapeMin));
					}
				}
				else if(mTimeBrushMode == SMOOTHEN)
				{
					//TODO: implement correct smoothing for boundary vertices
					Eigen::VectorXd smoothed_time_landscape = mSelectedSimHandler->time_landscape;
					Eigen::VectorXd count_edges = Eigen::VectorXd::Ones(smoothed_time_landscape.size());
					for (int i = 0; i < smp->F.rows(); ++i)
						for (int j = 0; j < 3; ++j)
						{
							smoothed_time_landscape(smp->F(i, j)) += (mSelectedSimHandler->time_landscape(smp->F(i, next3(j))) + mSelectedSimHandler->time_landscape(smp->F(i, prev3(j)))) / 2;
							count_edges(smp->F(i, j)) += 1;
						}
					smoothed_time_landscape.array() /= count_edges.array();
					for (int i = 0; i < smp->uv.rows(); ++i)
					{
						double dist = (picked_point.head(2).transpose() - smp->uv.row(i)).norm();
						if (mTimeLandscapeRadius < dist) continue;
						mSelectedSimHandler->time_landscape(i) = smoothed_time_landscape(i);
					}
				}
				mSelectedSimHandler->update_TimeLandscapeBatch();
			}

			u_last_rightbutton_update = getElapsedSeconds();
		}
	}

	// Internal update

	if (!mSimulationHandlers.empty())
	{
		static double u_time_last = getElapsedSeconds();
		if (mPlayback)
		{
			double time_passed_scaled = (getElapsedSeconds() - u_time_last) * mPlaySpeed;
			if (mSelectedSimHandler->simulator->get_timestep() <= time_passed_scaled)
			{
				int skip_frames = floor(time_passed_scaled / mSelectedSimHandler->simulator->get_timestep());
				mCurrentFrame += skip_frames;
				u_time_last += skip_frames / mPlaySpeed * mSelectedSimHandler->simulator->get_timestep();
				if (mPlayLoop && mSelectedSimHandler->simulator->get_simulated_frames() - 1 < mCurrentFrame) mCurrentFrame = (mLockRestConfiguration ? 1 : 0);
				mCurrentFrame = min(mCurrentFrame, mSelectedSimHandler->simulator->get_simulated_frames() - 1);
			}
			mCurrentFrameAlpha = time_passed_scaled / mSelectedSimHandler->simulator->get_timestep();
			mCurrentFrameAlpha = mCurrentFrameAlpha - floor(mCurrentFrameAlpha);
			if (mCurrentFrame == mSelectedSimHandler->simulator->get_simulated_frames() - 1) mCurrentFrameAlpha = 0;
			if (mCurrentFrame == 0) mCurrentFrameAlpha = 0; // no smoothing between rest configuration and simulation
		}
		else
		{
			u_time_last = getElapsedSeconds();
			mCurrentFrameAlpha = 0;
		}

		static int u_CurrentDisplayFrame = -1;
		if (!mSimulationHandlers.empty())
		{
			if (u_CurrentDisplayFrame != mCurrentFrame || mPlaySmooth)
			{
				u_CurrentDisplayFrame = mCurrentFrame;
				for (auto& shandler : mSimulationHandlers)
				{
					shandler->update_RigidBodyBatch(mCurrentFrame, mCurrentFrameAlpha);
					shandler->update_MembraneBatch(mCurrentFrame, mCurrentFrameAlpha);
				}
				update_ActuatorsBatch();
			}
		}
	}

	/*
	gl::VboMeshRef mVboMesh = mPrimitiveWireframe->getVboMesh();

	auto mappedPosAttrib = mVboMesh->mapAttrib3f(geom::Attrib::POSITION, false);
	for (int i = 0; i < mVboMesh->getNumVertices(); i++)
	{
		vec3 &pos = *mappedPosAttrib;
		mappedPosAttrib->x = psys->V(i,0);
		mappedPosAttrib->y = psys->V(i,1);
		mappedPosAttrib->z = psys->V(i,2);
		++mappedPosAttrib;
	}
	mappedPosAttrib.unmap();
	*/

	showGui();
}

void SmpApp::draw()
{
	// Prepare for drawing.
	gl::clear();
	gl::setMatrices(*mCameraActive);

	// Enable the depth buffer.
	gl::enableDepthRead();
	gl::enableDepthWrite();

	gl::enableFaceCulling();

	//gl::ScopedColor colorScope(Color(1, 1, 1));

	// Rotate it slowly around the y-axis.
	gl::ScopedModelMatrix matScope;
	if (mRotation) gl::rotate(float(getElapsedSeconds() / 5), 0, 0, 1);
	
	gl::ScopedTextureBind scopedTextureBind(mTexture);
	mPhongShader->uniform("uTexturingMode", mTexturingMode);
	if (mTexturingMode == PROCEDURAL) mPhongShader->uniform("uFreq", ivec2(1, 1));
	else if (mTexturingMode == SAMPLER) mPhongShader->uniform("uFreq", ivec2(5, 5));

	if (mShowCoordinateFrame) gl::drawCoordinateFrame(10);

	for (auto& shandler : mSimulationHandlers)
	{
		shandler->draw_solid();
		if (shandler->visible && mShowAlign)
		{
			auto& align = shandler->simulator->get_smp()->base_align;
			auto V = shandler->simulator->get_positions_at_frame(mCurrentFrame, mCurrentFrameAlpha);
			for (auto i : align)
			{
				float r = shandler->simulator->get_smp()->options.base_th + shandler->simulator->get_smp()->options.base_p + 0.3f;
				ColorA color = shandler->color;
				color.r = 1.0 - color.r;
				color.g = 1.0 - color.g;
				color.b = 1.0 - color.b;
				color.a = 0.6;
				gl::color(color);
				gl::cullFace(GL_BACK);
				gl::Batch::create(geom::Sphere(Sphere(vec3(V(i, 0), V(i, 1), V(i, 2)), r)), mSelectionShader)->draw();
			}
		}
	}
	if (mSelectedSimHandler && mShowTimeLandscapeEditor) mSelectedSimHandler->draw_time_landscape();
	if (mSelectedSimHandler && mSelectedSimHandler->visible && mActuatorsBatch)
	{
		mWireframeShader->uniform("uBrightness", 1.0f);
		mWireframeShader->uniform("uOpaqueness", 1.0f);
		gl::cullFace(GL_BACK);
		mActuatorsBatch->draw();
	}
	for (auto& shandler : mSimulationHandlers)
	{
		shandler->draw_stencil();
		shandler->draw_membrane();
	}

	for (auto& shandler : mSimulationHandlers)
	{
		if (!shandler->visible) continue;
		Eigen::MatrixXd a, b;
		Eigen::VectorXi k;
		if (mShowLinearSprings)
		{
			shandler->simulator->get_springs_at_frame(a, b, k, mCurrentFrame, mPlaySmooth ? mCurrentFrameAlpha : 0);
			for (int i = 0; i < a.rows(); ++i)
			{
				vector<vec3> edg;
				edg.push_back(vec3(a(i, 0), a(i, 1), a(i, 2)));
				edg.push_back(vec3(b(i, 0), b(i, 1), b(i, 2)));
				gl::color(1, 1 - k(i), 0);
				gl::draw(edg, false);
			}
		}
		if (mShowScissorSprings)
		{
			shandler->simulator->get_scissors_at_frame(a, b, mCurrentFrame, mPlaySmooth ? mCurrentFrameAlpha : 0);
			for (int i = 0; i < a.rows(); ++i)
			{
				vector<vec3> edg;
				edg.push_back(vec3(a(i, 0), a(i, 1), a(i, 2)));
				edg.push_back(vec3(b(i, 0), b(i, 1), b(i, 2)));
				gl::color(0, 0.6, 0);
				gl::draw(edg, false);
			}
		}
		if (mShowTimeLandscapeEditor && mShowTimeRanges && 0 < shandler->time_ranges.size())
		{
			const SmpData* smp = shandler->simulator->get_smp();
			for (int i = 0; i < smp->FF.rows(); ++i)
			{
				if (0 <= mPickedActuator && i != mPickedActuator) continue;
				for (int j = 0; j < 3; ++j)
				{
					vector<vec3> edg;
					Eigen::Vector2d v2d = 0.5 * (smp->Fc2d.row(smp->FF(i, 0)) + smp->Fc2d.row(smp->FF(i, 1)));
					edg.push_back(vec3(v2d(0) + j * 0.05 - 0.05, v2d(1), j < 2 ? shandler->time_ranges(i, 0 + j * 2) : 0));
					edg.push_back(vec3(v2d(0) + j * 0.05 - 0.05, v2d(1), j < 2 ? shandler->time_ranges(i, 1 + j * 2) : 120));
					if (j < 2) gl::color(1 - j * 0.5, 0.3 + j * 0.5, 1);
					else gl::color(0.5, 0.5, 0.5);
					gl::draw(edg, false);
				}
			}
		}
	}

#ifdef USE_GEOGRAM
	if (mShowSurfaceImport) mSmpRemeshHandler.draw();
#endif

	draw_grid();

	gl::disableDepthWrite();
	gl::disableDepthRead();
		
#ifdef USE_GEOGRAM
	mSmpRemeshHandler.draw_2d();
#endif
}

void SmpApp::resize()
{
	mCameraPersp.setAspectRatio(getWindowAspectRatio());
	mCameraOrtho.setOrtho(-mOrthoWidth, mOrthoWidth, -mOrthoWidth / getWindowAspectRatio(), mOrthoWidth / getWindowAspectRatio(), 0.1, 1000);

	if (mWireframeShader) mWireframeShader->uniform("uViewportSize", vec2(getWindowSize()));
}

void SmpApp::keyDown( KeyEvent event )
{
	if (event.isControlDown())
	{
		switch (event.getCode()) {
		case KeyEvent::KEY_o:
			read_file_info();
			model_open_stencil_request = true;
			break;
		}
	}
	else if (event.isAltDown())
	{
		switch (event.getCode()) {
		case KeyEvent::KEY_x:
			mExitRequest = true;
			break;
		}
	}
	else
	{
		switch (event.getCode()) {
		case KeyEvent::KEY_LEFT:
			if (mSelectedSimHandler)
			{
				mCurrentFrame -= 1;
				mCurrentFrame = max(mCurrentFrame, mLockRestConfiguration ? 1 : 0);
			}
			break;
		case KeyEvent::KEY_RIGHT:
			if (mSelectedSimHandler)
			{
				mCurrentFrame += 1;
				mCurrentFrame = min(mCurrentFrame, mSelectedSimHandler->simulator->get_simulated_frames() - 1);
			}
			break;
		case KeyEvent::KEY_SPACE:
			mPlayback = !mPlayback;
			break;
		}
	}
}


void SmpApp::mouseDown(MouseEvent event)
{
	mRecenterCamera = false;
	if (mCameraActive == &mCameraPersp) mCamUi.mouseDown(event);
	else if (mCameraActive == &mCameraOrtho)
	{
		if (event.isLeft()) mCameraOrthoMousePos = event.getPos();
	}
	if (event.isLeft())
	{
		mLastMouseLeftDownTime = getElapsedSeconds();
		mMouseLeftButtonDown = true;
		mMouseLeftClickDown = true;
	}
	if (event.isRight())
	{
		mLastMouseRightDownTime = getElapsedSeconds();
		mMouseRightButtonDown = true;
		mMouseRightClickDown = true;
	}
}


void SmpApp::mouseUp(MouseEvent event)
{
	if (mCameraActive == &mCameraPersp) mCamUi.mouseUp(event);
	if (getElapsedSeconds() - mLastMouseLeftDownTime < 0.2f && event.isLeft() && mMouseLeftClickDown)
	{
		if (!mSimulationHandlers.empty())
		{
			int picked_actuator_id = performActuatorPicking();
			int picked_base_id = -1;
			if (picked_actuator_id == -1 && mShowAlign)
			{
				picked_base_id = performBasePicking();
				if (-1 < picked_base_id)
				{
					auto& align = mSelectedSimHandler->simulator->get_align();
					if (align.erase(picked_base_id) == 0) align.insert(picked_base_id);
				}
			}
			if (picked_base_id == -1)
			{
				mPickedActuator = (picked_actuator_id == mPickedActuator ? -1 : picked_actuator_id);
				update_ActuatorsBatch();
			}
		}
	}
	if (event.isLeft())
	{
		mMouseLeftButtonDown = false;
		mMouseLeftClickDown = false;
	}
	if (event.isRight())
	{
		mMouseRightButtonDown = false;
		mMouseRightClickDown = false;
	}
}

void SmpApp::mouseDrag(MouseEvent event)
{
	if (mMousePos == event.getPos()) return;
	mMousePos = event.getPos();
	if (mCameraActive == &mCameraPersp && !event.isRightDown()) mCamUi.mouseDrag(event);
	if (mCameraActive == &mCameraOrtho && event.isLeftDown())
	{
		//TODO: only works for top view
		float u = (mMousePos.x - mCameraOrthoMousePos.x) / (float)getWindowWidth();
		float v = (mMousePos.y - mCameraOrthoMousePos.y) / (float)getWindowHeight();
		mCameraOrthoMousePos = mMousePos;
		vec3 eye = mCameraOrtho.getEyePoint();
		eye[0] -= u * mOrthoWidth * 2;
		eye[1] += v * mOrthoWidth * 2 / getWindowAspectRatio();
		mCameraOrtho.setEyePoint(eye);
	}
	mMouseLeftClickDown = false;
	mMouseRightClickDown = false;
}

void SmpApp::mouseMove(MouseEvent event)
{
	mMousePos = event.getPos();
}

void SmpApp::mouseWheel(MouseEvent event)
{
	if (mCameraActive == &mCameraPersp) mCamUi.mouseWheel(event);
	else
	{
		mOrthoWidth -= event.getWheelIncrement() * mOrthoWidth / 5;
		resize(); //TODO: dirty implementation
	}
}

template<typename DerivedA, typename DerivedB>
Eigen::Vector3d SmpApp::pick_triangle(const Eigen::DenseBase<DerivedA>& V, const Eigen::DenseBase<DerivedB>& F, int& picked)
{
	Eigen::Vector3d picked_point = Eigen::Vector3d::Zero();

	// Generate a ray from the camera into our world. Note that we have to
	// flip the vertical coordinate.
	float u = mMousePos.x / (float)getWindowWidth();
	float v = mMousePos.y / (float)getWindowHeight();
	Ray ray;
	if (mCameraActive == &mCameraPersp) ray = mCameraActive->generateRay(u, 1.0f - v, mCameraActive->getAspectRatio());
	else
	{
		//TODO: this is valid only for top view
		vec3 eye = mCameraOrtho.getEyePoint();
		eye[0] += mOrthoWidth * (2 * u - 1);
        eye[1] += mOrthoWidth * (1 - 2 * v) / getWindowAspectRatio(); //TODO: not ideal implementation of camera aspect ratio (cinder function fails)
		ray.setOrigin(eye);
		ray.setDirection(vec3(0, 0, -1));
	}

	//	// The coordinates of the bounding box are in object space, not world space,
	//	// so if the model was translated, rotated or scaled, the bounding box would not
	//	// reflect that. One solution would be to pass the transformation to the calcBoundingBox() function:
	//	AxisAlignedBox worldBoundsExact = mTriMesh->calcBoundingBox(); // slow
	//
	//	drawCube(mObjectBounds, Color(1, 1, 0));
	//
	//	// Draw the exact bounding box in orange.
	//	drawCube(worldBoundsExact, Color(1, 0.5f, 0));
	//
	//	// Draw the approximated bounding box in cyan.
	//	drawCube(worldBoundsApprox, Color(0, 1, 1));
	//
	//	// Perform fast detection first - test against the bounding box itself.
	//	if (!worldBoundsExact.intersects(ray))
	//		return false;

	// Set initial distance to something far, far away.
	float result = FLT_MAX;

	float distance = 0.0f;
	picked = -1;
	for (size_t i = 0; i < F.rows(); ++i)
	{
		vec3 v0(V(F(i, 0), 0), V(F(i, 0), 1), V(F(i, 0), 2));
		vec3 v1(V(F(i, 1), 0), V(F(i, 1), 1), V(F(i, 1), 2));
		vec3 v2(V(F(i, 2), 0), V(F(i, 2), 1), V(F(i, 2), 2));

		if (ray.calcTriangleIntersection(v0, v1, v2, &distance))
		{
			if (distance < result)
			{
				result = distance;
				picked = i;
			}
		}
	}

	// Did we have a hit?
	if (distance > 0) {
		vec3 picked_pnt = ray.calcPosition(result);
		picked_point(0) = picked_pnt.x;
		picked_point(1) = picked_pnt.y;
		picked_point(2) = picked_pnt.z;

		// Barycentric
		//Eigen::Vector3d r{ picked_pnt[0], picked_pnt[1], picked_pnt[2] };
		//Eigen::MatrixXd T(3, 2);
		//T.col(0) = V.row(F(picked, 0)) - V.row(F(picked, 2));
		//T.col(1) = V.row(F(picked, 1)) - V.row(F(picked, 2));
		//Eigen::Vector2d lam = T.householderQr().solve(r - V.row(F(picked, 2)).head(3).transpose());
		//bary.head(2) = lam;
		//bary(2) = 1 - lam.sum();
	}
	return picked_point;
}

int SmpApp::performActuatorPicking()
{
	if (mSimulationHandlers.empty()) return -1;
	if (!mSelectedSimHandler->visible) return -1;
	
	auto V = mSelectedSimHandler->simulator->get_positions_at_frame(mCurrentFrame); // not accounting smooth playback
	auto F = mSelectedSimHandler->simulator->get_actuator_faces();

	int picked_face;
	pick_triangle(V, F, picked_face);
	if (picked_face < 0)
	{
		mActuatorsBatch.reset();
		return -1;
	}

	update_ActuatorsBatch();
		
	return F(picked_face, 3);
}

int SmpApp::performBasePicking()
{
	if (mSimulationHandlers.empty()) return -1;
	if (!mSelectedSimHandler->visible) return -1;

	auto  V = mSelectedSimHandler->simulator->get_body_vertices_at_frame(mCurrentFrame); // not accounting smooth playback
	auto& F = mSelectedSimHandler->simulator->get_body_faces_indexed();

	int picked_face;
	pick_triangle(V, F, picked_face);
	if (picked_face < 0) return -1;

//	update_ActuatorsBatch();

	return F(picked_face, 3);
}

void SmpApp::fileDrop( FileDropEvent event )
{
	try {
		fs::path path(event.getFile(0));
		read_file_info(path);
		model_open_stencil_request = true;
	}
	catch( const exception &exc ) {
	}
}

void SmpApp::update_ActuatorsBatch()
{
	if (mPickedActuator < 0)
	{
		mActuatorsBatch.reset();
		return;
	}

	Eigen::MatrixXd V = mSelectedSimHandler->simulator->get_positions_at_frame(mCurrentFrame, mCurrentFrameAlpha);

	TriMesh::Format fmt = TriMesh::Format().positions().colors(4);
	TriMeshRef mActuatorsTriMesh = TriMesh::create(fmt);

	Eigen::MatrixXi Fa = mSelectedSimHandler->simulator->get_actuator_faces();
	for (int i = 0; i < V.rows(); ++i)
	{
		glm::vec3 pos(V(i, 0), V(i, 1), V(i, 2));
		mActuatorsTriMesh->appendPosition(pos);
		mActuatorsTriMesh->appendColorRgba(ColorA(0.2, 0.2, 1, 0.5));
	}
	for (int i = 0; i < Fa.rows(); ++i)
	{
		if (Fa(i, 3) == mPickedActuator)
			mActuatorsTriMesh->appendTriangle(Fa(i, 0), Fa(i, 1), Fa(i, 2));
	}

	mActuatorsBatch = gl::Batch::create(*mActuatorsTriMesh.get(), mSelectionShader);
}

CINDER_APP( SmpApp, RendererGl( RendererGl::Options().msaa( 16 ) ), prepareSettings )
