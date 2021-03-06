#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include "SmartPtr.h"

auto_ptr<MyShareable> gime_mine() { return mine; }
auto_ptr<MyShareable>* gime_mine_ptr() { return &mine; }
auto_ptr<MyShareable>& gime_mine_ref() { return mine; }

void pass_mine_sp(auto_ptr<MyShareable> )
//{ cout << "Use count " << p->use_count() << "\n"; }
{ } // cout << "pass_mine_sp: underlying ptr " << p.get() << "\n"; }

void pass_mine_sp_ref(auto_ptr<MyShareable>& )
//{ cout << "Use count " << p->use_count() << "\n"; }
{ } // cout << "pass_mine_sp_ref: underlying ptr " << p.get() << "\n"; }

void pass_mine_sp_ptr(auto_ptr<MyShareable>* )
//{ cout << "Use count " << p->use_count() << "\n"; }
{ } // cout << "pass_mine_sp_ptr: underlying ptr " << p->get() << "\n"; }

void pass_mine_rp(MyShareable) {}
void pass_mine_rp_ref(const MyShareable&) {}
void pass_mine_rp_ptr(const MyShareable*) {}

#pragma GCC diagnostic pop
