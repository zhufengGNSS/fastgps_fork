/*
* Copyright (c) 2008, Scott Gleason and Morgan Quigley
* All rights reserved.
*
* Written by Scott Gleason and Morgan Quigley
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the authors' names nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "wx/wx.h" 
#include "fastgps_frame.h"
#include "fastgps.h"

static FastGpsFrame *g_main_frame;
static FILE         *g_message_log;

enum
{
  ID_Quit = 1,
  ID_About,
  ID_FASTGPS_START_BUTTON,
  ID_FASTGPS_LOG_BOX,
};

BEGIN_EVENT_TABLE(FastGpsFrame, wxFrame)
  EVT_MENU(ID_Quit, FastGpsFrame::OnQuit)
  EVT_MENU(ID_About, FastGpsFrame::OnAbout)
  EVT_BUTTON(ID_FASTGPS_START_BUTTON, FastGpsFrame::startButtonClick)
END_EVENT_TABLE() 

FastGpsFrame::FastGpsFrame(const wxString& title, const wxPoint& pos, 
                           const wxSize& size)
: wxFrame((wxFrame *)NULL, -1, title, pos, size)
{
  g_main_frame = this;
	wxBoxSizer *topsizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *vertsizer = new wxBoxSizer(wxVERTICAL);
	wxStaticText *text = new wxStaticText(this, -1, _T("fastgps"));
	topsizer->Add(vertsizer, 1, wxALIGN_CENTER | wxALL, 10);
	vertsizer->Add(text, 1, wxALIGN_CENTER | wxALL, 10);
	SetSizer(topsizer);
  wxMenu *menuFile = new wxMenu;
  menuFile->Append( ID_About, _T("&About...") );
  menuFile->AppendSeparator();
  menuFile->Append( ID_Quit, _T("E&xit") ); 
  wxMenuBar *menuBar = new wxMenuBar;
  menuBar->Append( menuFile, _T("&File") ); 
  SetMenuBar( menuBar );
  CreateStatusBar();
  SetStatusText( _T("Welcome to the fastgps software receiver!") );  
  //     Button1 = new wxButton(this, ID_BUTTON1, _T("Start"), wxPoint(240,232), wxDefaultSize, 0, wxDefaultValidator, _T("ID_BUTTON1"));
  start_button = new wxButton(this, ID_FASTGPS_START_BUTTON, _T("Start"), 
                              wxPoint(40,125));
  log_box = new wxTextCtrl(this, ID_FASTGPS_LOG_BOX, wxEmptyString, 
                           wxPoint(140,24), wxSize(456,292), wxTE_MULTILINE, 
                           wxDefaultValidator, _T("ID_FASTGPS_LOG_BOX"));
}
 
void FastGpsFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
  Close(TRUE);
} 

void FastGpsFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
  wxMessageBox(_T("This is the fastgps software receiver"),
      _T("About fastgps"), wxOK | wxICON_INFORMATION, this);
}

void FastGpsFrame::startButtonClick(wxCommandEvent& WXUNUSED(event))
{
  log_box->WriteText(_T("fastgps stating"));
  g_message_log = fopen("fastgps_messages.txt","w");
  run_fastgps();
  fclose(g_message_log);
}
/*
void FastGpsFrame::AddMessage(char *msg)
{
  wxString bbb(msg, wxConvUTF8);
  log_box->WriteText(bbb);  // make sure  wxTextCtrl is set to MULTILINE=true (expand EditStyle)
  wxYieldIfNeeded();
  fprintf(message_log, msg); // log message to file
}
*/
void fastgps_printf(const char *format, ...)
{
  va_list args;
  va_start(args, format);
  char msgbuf[65535]; // that's one long printf
  vsprintf(msgbuf, format, args);
  va_end(args);
  wxString wx_str(msgbuf, wxConvUTF8);
  g_main_frame->log_box->WriteText(wx_str); 
  wxYieldIfNeeded();
  fprintf(g_message_log, msg); // log message to file
}

