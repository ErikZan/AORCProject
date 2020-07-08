ground := Matrix(4, 4, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], attributes =
[protected, protected, ground]);
omegaMV := table([(frame)=(Matrix(4, 4, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1
]], attributes = [protected, protected, ground])),(obj)=VECTOR,(comps)=(Matrix(\
4, 1, [[0],[0],[0],[0]]))]);
draw_link := proc (P1::POINT, P2::POINT, data, sol, wb := .17e-1, {col := 
"SkyBlue"}, $) local tmp1, tmp2; tmp1 := simplify(subs(sol,data,[
MBSymba_r6_kinematics:-comp_XYZ(P1,ground)])); tmp2 := simplify(subs(sol,data,[
MBSymba_r6_kinematics:-comp_XYZ(P2,ground)])); plottools[plottools:-arrow](tmp1
,tmp2,wb,0,0,cylindrical_arrow,color = col); end proc;
small_vars := {};
dependent_coordinates := false;
draw_quadcopter := proc (data, sol) draw_sphere(Obody,r = 
MBSymba_r6_foundations:-`*`(10,.5e-1),data,sol,col = "LightSalmon"), 
draw_sphere(OmF,r = MBSymba_r6_foundations:-`*`(10,.2e-1),data,sol), 
draw_sphere(OmB,r = MBSymba_r6_foundations:-`*`(10,.2e-1),data,sol), 
draw_sphere(OmL,r = MBSymba_r6_foundations:-`*`(10,.2e-1),data,sol), 
draw_sphere(OmR,r = MBSymba_r6_foundations:-`*`(10,.2e-1),data,sol), draw_link(
OmF,OmB,data,sol,col = "Blue"), draw_link(OmL,OmR,data,sol,col = "Blue"); end 
proc;
fname := "Functions.mpl";
mframe_flag := false;
transform_obj := proc (RF::frame, obj::function, col::string := "", transp::
scalar := 0, $) local FF, ff, f, ip, obj_new; if col <> "" then ip := ListTools
[Search](COLOUR,map2(op,0,[op(obj)])); obj_new := subsop(ip = COLOUR(RGB,op(
ColorTools[ToRGB24](col))),obj); else obj_new := obj; end if; ip := ListTools[
Search](TRANSPARENCY,map2(op,0,[op(obj_new)])); obj_new := subsop(ip = 
TRANSPARENCY(transp),obj_new); FF := evalf(simplify(RF) . <x, y, z, 1>); 
convert(FF[1 .. 3],list); ff := unapply(%,[x, y, z]); f := plottools[plottools
:-transform]((x, y, z) -> ff(x,y,z)); f(obj_new); end proc;
master_frame := Matrix(4, 4, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 
attributes = [protected, protected, ground]);
draw_sphere := proc (P::POINT, data, sol, {col := "Goldenrod", r := .1e-1}, $)
local tmp; tmp := subs(sol,data,[MBSymba_r6_kinematics:-comp_XYZ(P,ground)]); 
plottools[plottools:-sphere](tmp,r,color = col); end proc;
draw_RF := proc (RF::frame, data, sol, col, scale, $) local p, vX, vY, vZ, RFd;
RFd := subs(sol,data,RF); p := [MBSymba_r6_kinematics:-comp_XYZ(
MBSymba_r6_kinematics:-origin(RFd),ground)]; vX := [MBSymba_r6_kinematics:-
comp_XYZ(MBSymba_r6_kinematics:-uvec_X(RFd),ground)]; vY := [
MBSymba_r6_kinematics:-comp_XYZ(MBSymba_r6_kinematics:-uvec_Y(RFd),ground)]; vZ
:= [MBSymba_r6_kinematics:-comp_XYZ(MBSymba_r6_kinematics:-uvec_Z(RFd),ground)]
; plots[plottools:-arrow](p,MBSymba_r6_foundations:-`*`(vX,scale),shape = 
cylindrical_arrow,color = "Red",fringe = col), plots[plottools:-arrow](p,
MBSymba_r6_foundations:-`*`(vY,scale),shape = cylindrical_arrow,color = "Green"
,fringe = col), plots[plottools:-arrow](p,MBSymba_r6_foundations:-`*`(vZ,scale)
,shape = cylindrical_arrow,color = "Blue",fringe = col); end proc;
