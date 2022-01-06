# Linear Systems and Matrices

## Gaussian Elimination

The technique studied here makes solving large systems of equations (e.g. 200 equations in 200 unknowns) practical, and it can be implemented on a computer in a basic programming language. You might be surprised to learn that it is necessary to solve vast systems of linear equations as a matter of routine in many practical scientific applications.

### Motivation

Gaussian elimination is a systematic technique for solving systems of linear equations, which are of the form

\begin{align}
a_{1,1} x_1 + a_{1,2} x_2 + a_{1,3} x_3 + \dots + a_{1,n} x_n &= b_1 \\

   a_{2,1} x_1 + a_{2,2} x_2 + a_{2,3} x_3 + \dots + a_{2,n} x_n &= b_2\\
                                                                     \vdots \\
 a_{m,1} x_1 + a_{m,2} x_2 + a_{m,3} x_3 + \dots + a_{m,n} x_n &= b_m
\end{align}


where $a_{i,j}$ are constants. Usually there are the same number of equations as unknowns, so $m = n$.  If $m < n$ then the system is undetermined, if $m > n$, the system is over constrained.

```{admonition} Example
:class: tip
\begin{align}4 x_1 - x_2 &= 1\\ -2 x_1 + 3 x_2 &= 12 \end{align}
Each equation here defines a line, and we are looking for a point which satisfies both equations, which means that the lines intersect.

From the first line we obtain $x_2 = 4 x_1 - 1$, and by substituting this into the second line, we obtain $x = \frac{3}{2}$, $y = 5$. Two equations with two unknowns will always give a unique solution, unless the lines are parallel (and so the equations are just a scaling of each other).

-	If they are parallel and distinct, there are no solutions because there are no points that lie on both lines.
-	If they are parallel and coincident (same line), there are an infinite number of solutions.
```
3.1.2. A systematic technique for solving systems of equations
We will begin by finding a solution to the following system:

The equations have been labelled $r_1,r_2,r_3$
First, we will use $r_1$ to eliminate $x_1$ from $r_2$ and $r_3$. This gives two equations in two unknowns. Then, we will use $r_2$ to eliminate $x_2$ from $r_3$.
The steps are written out below:    

The solution for $x_3$ can now be read off from $r_3$, $x_2$ can be obtained from $r_2$ using the result for $x_3$ and $x_1$ can be obtained from $r_1$ using the results for $x_1$ and $x_2$ - This is known as “back-substitution”.  

These manipulations can be conveniently done by looking only at the coefficients, which we collect together in a form called the augmented matrix:

We can see that the algorithm (described in the box below) works by eliminating the coefficients below the leading diagonal, which is highlighted in pink in the image.

Box 3.1 Naive Gaussian elimination algorithm (obtaining upper triangular form)
(A). Choose initial pivot
We choose the first element from the leading diagonal as the pivot element
(B). Row reduction step
Add multiples of the pivot row to the rows below, to obtain zeros in the pivot column below the leading diagonal.
(C). Choose new pivot
The pivot moves to the next element on the leading diagonal
Repeat from step (B) until the matrix is in upper triangular form (containing all zeros below the leading diagonal).
The solutions can then be obtained by back-substitution

Examples for you to try:
You can use the following app to generate a system of three equations in three unknowns and solve it using the Gaussian elimination algorithm given in Box 3.1. In the solutions, the pivot element in each step is highlighted in blue.
https://www.wolframcloud.com/obj/ucqssjm/Published/gaussian-elimination


3.1.3. Generalisation
The naive algorithm introduced here can be generalised to include additional row operations. In general, the acceptable row operations that we can perform are:
•	multiplication of any row by a constant
•	addition of (a multiple of) any row to any other
•	swapping any two rows
It is often possible to apply these steps creatively to get a result with greater efficiency than using the naive algorithm described above. However, in this course the assessed questions will all require you to carry out the naive algorithm. There are two reasons for this:
1.	To test your ability to follow an algorithm or set of rules. Once you can follow an algorithm precisely, accurately and definitively, this is the starting point at which you could programme a computer to automate the method.

2.	For ease of process-checking - for example, when using automated quiz-marking, or when comparing solutions.

It is also not necessary to stop at upper triangular form. Once the last row has been fully simplified, it can be used to obtain zeros above the main diagonal in the last column. Then, the second-last row is used to obtain zeros in the second-last column above the main diagonal, and so-on until the only non-zero elements remaining are on the main diagonal. Then the solutions can be simply read off from each row.
For instance, continuing with the naive row reduction for the example shown in 3.1.2, we obtain

We have obtained row-reduced form and the solutions for $x_1,x_2,x_3$ can now be read off from the final column.

Advanced application: Kirchoff's Law
In general, Gaussian elimination can be used to solve the problems obtained by applying Kirchoff's laws. For example, see the problems given here:
https://www.intmath.com/matrices-determinants/6-matrices-linear-equations.php
In the example below, the technique is applied to a case where the system of equations is under-determined, so a unique solution cannot be obtained.

Kirchoff's law states that current in = current out.
For the system of 4 nodes shown above, this gives us four equations:
\begin{align}y_3&=y_1+y_4\\y_1&=y_2+y_5\\y_2&=y_3+y_6\\y_4+y_5+y_6&=0\end{align}
The equations at each node are of the form $$c_1 y_1 +c_2 y_2+\dots +c_6 y_6=d$$. Written in augmented matrix form, the system is:
$$\left(\begin{array}{cccccc:c}1 & 0 & -1 & 1 & 0 & 0 &0\\-1 & 1 & 0 & 0 & 1 & 0&0 \\0 & -1 & 1 & 0 & 0 & 1&0 \\0 & 0 & 0 & -1 & -1 & -1 &0\\\end{array}\right)$$
The first 6 columns are for the coefficients of $$y_1,\dots y_6$$ and the last column is the constant term that appears on the right hand side of the equation.
Here, Gaussian elimination is applied as usual, but the complicating factor is that there are more unknowns than equations and so the system is under-determined. It has "free variables", that can be made to take any value!
After applying Gaussian elimination to the pivot elements in the first and second column, we can't do anything with the third column without spoiling the progress we've made in the first two columns, so we leave that one and move on to use the fourth column as a pivot. After that, we can't make any more progress so we stop. We obtain:
$$\left(\begin{array}{cccccc:c}1 & 0 & -1 & 0 & -1 & -1 & 0 \\0 & 1 & -1 & 0 & 0 & -1 & 0 \\0 & 0 & 0 & 1 & 1 & 1 & 0 \\0 & 0 & 0 & 0 & 0 & 0 & 0 \\\end{array}\right)$$
The first, second, and fourth columns are the pivot columns, and the other columns are all "free" since they can be obtained from a combination of the other columns.
Choosing $$(y_3,y_5,y_6)=(1,0,0)$$ gives the special solution $$s_1:(1,1,1,0,0,0)$$
Choosing $$(y_3,y_5,y_6)=(0,1,0)$$ gives the special solution $$s_2:(1,0,0,-1,1,0)$$
Choosing $$(y_3,y_5,y_6)=(0,0,1)$$ gives the special solution $$s_3:(1,1,0,-1,0,1)$$
The full solution space consists of all possible linear combinations $$a s_1 + b s_2 +c s_3$$.
