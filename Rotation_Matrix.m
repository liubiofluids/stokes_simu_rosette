function m = Rotation_Matrix(theta, indx_axis)
    m=zeros(3,3);
    m(1,:)=[1, 0, 0];
    m(2,:)=[0, cos(theta), sin(theta)];
    m(3,:)=[0, -sin(theta), cos(theta)];
    m=shftmat(m, 1-indx_axis, 1-indx_axis);
end

