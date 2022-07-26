/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright (c) Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

private:

//
// Members for landmarks and object display and manipulation
//

/// Widget for display control of landmarks
Fl_Button *viewLandmarks;

Fl_Button *refineTags;

Fl_Check_Button *viewTagGrid;

/// Widgets for landmark browsing
Fl_Multi_Browser *targetLandmarkBrowser;
Fl_Multi_Browser *sourceLandmarkBrowser;

/// Widget for ROI
Fl_Button *viewROI;

/// Callbacks for landmarks and objects
static void cb_viewROI(Fl_Button*, void*);
static void cb_trackTAG(Fl_Button*, void*);
static void cb_viewTagGrid(Fl_Button*, void*);
#ifdef HAVE_VTK
Fl_Check_Button *viewObjectMovie;
Fl_Button *warpObject;
static void cb_viewObjectMovie(Fl_Button*, void*);
static void cb_loadObject(Fl_Button*, void*);
static void cb_warpObject(Fl_Button*, void*);
static void cb_removeObject(Fl_Button*, void*);
#endif
static void cb_addLandmark(Fl_Button*, void*);
static void cb_deleteLandmark(Fl_Button*, void*);
static void cb_toggleLandmark(Fl_Input*, void*);
static void cb_insertLandmark(Fl_Button*, void*);
static void cb_replaceLandmark(Fl_Button*, void*);
static void cb_editLandmark(Fl_Button*, void*);
static void cb_browseLandmark(Fl_Browser*, void*);
static void cb_viewLandmarks(Fl_Button*, void*);
static void cb_fitLandmarks(Fl_Button*, void*);
static void cb_loadTargetLandmarks(Fl_Button*, void*);
static void cb_loadSourceLandmarks(Fl_Button*, void*);
static void cb_saveTargetLandmarks(Fl_Button*, void*);
static void cb_saveSourceLandmarks(Fl_Button*, void*);

public:

/// Show the object control window
void ShowObjectControlWindow();

/// Update the object control window
void UpdateObjectControlWindow();

/// Initialize the object control window
void InitializeObjectControlWindow();
