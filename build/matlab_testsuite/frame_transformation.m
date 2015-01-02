function T = frame_transformation(DOF,theta)

if ~isvector(theta) 
    error('myApp:argChk', 'Second argument must be a vector')
end

A_link = [0;0;0.045;-0.045;0;0;0];
Alpha_link = [-pi/2;pi/2;-pi/2;pi/2;-pi/2;pi/2;0];
D_link = [0;0;0.55;0;0.3;0;0.06];

    function Tt = one_transf(A_f,alpha_f,D_f,Theta_f)
       % Tt  = zeros(4);
        Tt = [cos(Theta_f), -sin(Theta_f)*cos(alpha_f), sin(Theta_f)*sin(alpha_f), A_f*cos(Theta_f);
     sin(Theta_f), cos(Theta_f)*cos(alpha_f),-cos(Theta_f)*sin(alpha_f), A_f*sin(Theta_f);
     0,            sin(alpha_f),              cos(alpha_f),              D_f;
     0,            0,                         0,                         1];
    end

T = eye(4);

for i=1:DOF
 T = T*one_transf(A_link(i),Alpha_link(i),D_link(i),theta(i));
end

end
