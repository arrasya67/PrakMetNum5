% 2x1 + x2 - x3 = 3
% 4x1 + 3x2 + x3 = 9
% -2x1 + x2 + 2x3 = 4

clc;
clear;

% Matriks A dan vektor b
A = [2 1 -1; 4 3 1; -2 1 2];
b = [3; 9; 4];
n = length(b);

disp('SISTEM PERSAMAAN LINEAR');
disp('Matriks A:');
disp(A);
disp('Vektor b:');
disp(b);

% 1. ELIMINASI GAUSS

disp('1. ELIMINASI GAUSS');

A1 = A;
b1 = b;

% Forward Elimination
for k = 1:n-1
    for i = k+1:n
        m = A1(i,k)/A1(k,k);
        A1(i,k:n) = A1(i,k:n) - m*A1(k,k:n);
        b1(i) = b1(i) - m*b1(k);
    end
end

disp('Matriks segitiga atas [U|y]:');
disp([A1 b1]);

% Backward Substitution
x_gauss = zeros(n,1);
x_gauss(n) = b1(n)/A1(n,n);
for i = n-1:-1:1
    x_gauss(i) = (b1(i) - A1(i,i+1:n)*x_gauss(i+1:n))/A1(i,i);
end

disp('Solusi Eliminasi Gauss:');
fprintf('x1 = %.4f\n', x_gauss(1));
fprintf('x2 = %.4f\n', x_gauss(2));
fprintf('x3 = %.4f\n', x_gauss(3));
disp('');

% 2. ELIMINASI GAUSS-JORDAN
 disp('2. ELIMINASI GAUSS-JORDAN');
Aug = [A b];

% Gauss-Jordan Elimination
for k = 1:n
    % Normalisasi pivot
    Aug(k,:) = Aug(k,:)/Aug(k,k);
    % Eliminasi kolom k pada baris lain
    for i = 1:n
        if i != k
            m = Aug(i,k);
            Aug(i,:) = Aug(i,:) - m*Aug(k,:);
        end
    end
end

disp('Matriks identitas [I|x]:');
disp(Aug);

x_jordan = Aug(:,n+1);

disp('Solusi Gauss-Jordan:');
fprintf('x1 = %.4f\n', x_jordan(1));
fprintf('x2 = %.4f\n', x_jordan(2));
fprintf('x3 = %.4f\n', x_jordan(3));
disp('');


% 3. DEKOMPOSISI LU

disp('3. DEKOMPOSISI LU');

L = eye(n);
U = A;
% Forward Elimination untuk membentuk L dan U
for k = 1:n-1
    for i = k+1:n
        m = U(i,k)/U(k,k);
        L(i,k) = m;
        U(i,k:n) = U(i,k:n) - m*U(k,k:n);
    end
end

disp('Matriks L (segitiga bawah):');
disp(L);
disp('Matriks U (segitiga atas):');
disp(U);

disp('Verifikasi A = LU:');
disp(L*U);

% Forward Substitution: Ly = b
y = zeros(n,1);
y(1) = b(1)/L(1,1);
for i = 2:n
    y(i) = (b(i) - L(i,1:i-1)*y(1:i-1))/L(i,i);
end

disp('Vektor y (Ly = b):');
disp(y);

% Backward Substitution: Ux = y
x_lu = zeros(n,1);
x_lu(n) = y(n)/U(n,n);
for i = n-1:-1:1
    x_lu(i) = (y(i) - U(i,i+1:n)*x_lu(i+1:n))/U(i,i);
end

disp('Solusi Dekomposisi LU:');
fprintf('x1 = %.4f\n', x_lu(1));
fprintf('x2 = %.4f\n', x_lu(2));
fprintf('x3 = %.4f\n', x_lu(3));
disp('');


% PERBANDINGAN HASIL

disp('============================================================');
disp('PERBANDINGAN HASIL KETIGA METODE');
disp('============================================================');
disp('Metode          |    x1    |    x2    |    x3    |');
disp('-------------------------------------------------------');
fprintf('Gauss           | %8.4f | %8.4f | %8.4f |\n', x_gauss(1), x_gauss(2), x_gauss(3));
fprintf('Gauss-Jordan    | %8.4f | %8.4f | %8.4f |\n', x_jordan(1), x_jordan(2), x_jordan(3));
fprintf('Dekomposisi LU  | %8.4f | %8.4f | %8.4f |\n', x_lu(1), x_lu(2), x_lu(3));
disp('-------------------------------------------------------');


