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
// Members for registration display and manipulation
//

Fl_Value_Slider *targetBlurring;
Fl_Value_Slider *sourceBlurring;
Fl_Value_Slider *targetResolution;
Fl_Value_Slider *sourceResolution;
Fl_Value_Slider *targetPadding;
Fl_Counter *numberOfIterations;
Fl_Counter *numberOfBins;
Fl_Counter *numberOfSteps;
Fl_Counter *lengthOfSteps;
Fl_Choice  *similarityMeasure;
Fl_Check_Button *level[5];
Fl_Return_Button *runRegistration;

/// Similarity measure menu
static Fl_Menu_Item menu_similarityMeasure[];

/// Callback for reading paramter
static void cb_loadParameter(Fl_Menu_Bar*, void*);

/// Callback for reading paramter
static void cb_saveParameter(Fl_Menu_Bar*, void*);

/// Callback for parameter updates
static void cb_similarityMeasure(Fl_Widget *, char *) ;

/// Callback for parameter updates
static void cb_numberOfBins(Fl_Counter*, void*);

/// Callback for parameter updates
static void cb_numberOfIterations(Fl_Counter*, void*);

/// Callback for parameter updates
static void cb_numberOfSteps(Fl_Counter*, void*);

/// Callback for parameter updates
static void cb_lengthOfSteps(Fl_Counter*, void*);

// Callbacks for parameter updates
static void cb_numberOfLevels(Fl_Button *, char *);

// Callbacks for parameter updates
static void cb_targetBlurring(Fl_Value_Slider *, void *);

// Callbacks for parameter updates
static void cb_targetResolution(Fl_Value_Slider *, void *);

/// Callback for registration target intensity padding
static void cb_targetPadding(Fl_Value_Slider *, void *);

// Callbacks for parameter updates
static void cb_sourceBlurring(Fl_Value_Slider *, void *);

// Callbacks for parameter updates
static void cb_sourceResolution(Fl_Value_Slider *, void *);

// Callbacks for registrations
static void cb_startRegistration(Fl_Button*, void*);

// Callbacks for registrations
static void cb_stopRegistration(Fl_Button*, void*);

// Callbacks for registrations
static void cb_guessRegistrationParameter(Fl_Button*, void*);

public:

/// Update the registration control window
void UpdateRegistrationControlWindow();

/// Initialize the registration control window
void InitializeRegistrationControlWindow();

