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

#include <mirtk/Image.h>
#include <mirtk/Transformation.h>
#include <mirtk/Registration.h>

#include "Fl_RView.h"
#include "Fl_RViewUI.h"

extern Fl_RViewUI *rviewUI;

Fl_RView::Fl_RView(int x, int y, int w, int h, const char *name) : Fl_Gl_Window(x, y, w, h, name)
{
  v = new RView(w, h);
  // See https://www.fltk.org/doc-1.3/osissues.html#osissues_macos, section "OpenGL and 'retina' displays"
  #if FL_API_VERSION >= 10304
    Fl::use_high_res_GL(1);
  #endif
}

void Fl_RView::draw()
{
  if (!valid()) {
    v->Resize(pixel_w(), pixel_h());
  }
  v->Draw();
}

int Fl_RView::handle(int event)
{
  char buffer1[256], buffer2[256], buffer3[256], buffer4[256], buffer5[256];

  const int event_x = pixel_event_x();
  const int event_y = pixel_event_y();
  const int event_dy = pixel_event_dy();

  switch (event) {
  case FL_KEYBOARD:
  case FL_SHORTCUT:
    if (Fl::event_key() == FL_Delete) {
      v->ClearContour();
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if ((Fl::event_key() == FL_BackSpace) && (Fl::event_shift() != 0)) {
      v->ClearContour();
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if (Fl::event_key() == FL_BackSpace) {
      v->UndoContour();
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if (Fl::event_key() == FL_Enter) {
#ifdef HAS_SEGMENTATION_PANEL
      if (rviewUI->_selected == 0) {
        fl_alert("Please, select a label!\n");
        return 1;
      } else {
        v->FillContour(rviewUI->_id, 0);
        v->ClearContour();
        v->SegmentationLabelsOn();
        v->SegmentationUpdateOn();
        v->Update();
#ifdef HAS_VISUALISATION_PANEL
        rviewUI->UpdateImageControlWindow();
#endif
        rviewUI->update();
        this->redraw();
        return 1;
      }
#endif
    }
    if (Fl::event_key() == FL_Tab) {
#ifdef HAS_SEGMENTATION_PANEL
      v->RegionGrowContour(event_x, event_y);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
#endif
    }
    v->cb_keyboard(Fl::event_text()[0]);
    rviewUI->update();
    this->redraw();
    return 1;
    break;
  case FL_PUSH:
    if ((Fl::event_button() == 1) && (Fl::event_shift() != 0)) {
      v->AddContour(event_x, event_y, FirstPoint);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
#ifdef HAS_SEGMENTATION_PANEL
    if ((Fl::event_button() == 3) && (Fl::event_shift() != 0)) {
      v->FillArea(event_x, event_y);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
#endif
    if ((Fl::event_button() == 1) && (Fl::event_ctrl() != 0)) {
      v->UpdateROI1(event_x, event_y);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if ((Fl::event_button() == 3) && (Fl::event_ctrl() != 0)) {
      v->UpdateROI2(event_x, event_y);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if (Fl::event_button() == 1) {
      v->SetOrigin(event_x, event_y);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    break;
  case FL_DRAG:
    if ((Fl::event_state() & (FL_BUTTON1 | FL_SHIFT)) == (FL_BUTTON1 | FL_SHIFT)) {
      v->AddContour(event_x, event_y, NewPoint);
      v->Update();
      rviewUI->info_voxel->value(buffer1);
      rviewUI->info_world->value(buffer2);
      rviewUI->info_target->value(buffer3);
      rviewUI->info_source->value(buffer4);
      rviewUI->info_segmentation->value(buffer5);
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if ((Fl::event_state() & (FL_BUTTON1 | FL_CTRL)) == (FL_BUTTON1 | FL_CTRL)) {
      v->UpdateROI1(event_x, event_y);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if ((Fl::event_state() & (FL_BUTTON3 | FL_CTRL)) == (FL_BUTTON3 | FL_CTRL)) {
      v->UpdateROI2(event_x, event_y);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    break;
  case FL_RELEASE:
    if ((Fl::event_button() == 1) && (Fl::event_shift() != 0)) {
      v->AddContour(event_x, event_y, LastPoint);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    return 0;
  case FL_MOVE:
    v->MousePosition(event_x, event_y);
    rviewUI->update();
    this->redraw();
    return 1;
    break;
  case FL_MOUSEWHEEL:
    v->MouseWheel(event_x, event_y, event_dy);
    v->Update();
    rviewUI->update();
    this->redraw();
    return 1;
    break;
  case FL_FOCUS:
    return 1;
    break;
  case FL_UNFOCUS:
    return 1;
    break;
  case FL_ENTER:
    Fl::focus(this);
    return 1;
    break;
  case FL_LEAVE:
    return 1;
    break;
  default:
    return Fl_Gl_Window::handle(event);
  }
  return 0;
}

