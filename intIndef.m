function [ f1 ] = intIndef( f0, x )
%INTINDEF Summary of this function goes here
%   Detailed explanation goes here

f1 = @(f0,x,dx,varargin) sum(f0( 0:dx:x, varargin)*dx);

end

