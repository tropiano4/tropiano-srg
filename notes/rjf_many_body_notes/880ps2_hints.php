<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"> 

<HTML>

<HEAD>
  <TITLE>OSU Physics: Physics 880.05 Assignment #2 Hints</TITLE>
  <LINK REV="made" HREF="mailto:webmaster@www.physics.ohio-state.edu">
  <LINK REL=STYLESHEET HREF="26x_style.css" TYPE="text/css">  
</HEAD>

<BODY BGCOLOR="#FFFFFF">


<?
   include ("physics/common_short2.php");
   physics_header("","");
?>
            

<H1>Physics 880.05: Assignment #2</H1>

Here are some hints, suggestions, and comments on the assignment.

<OL>

 <LI>MATLAB Sandbox, Part II.
   <ol TYPE="a">
     <li>This is really just an informational problem.  You can
     learn more through the MATLAB help by looking up 
     <tt>rand</tt> and <tt>randn</tt>.  (Note: MATLAB help is
     generally really good!)
     I've learned that the latest MATLAB versions do the seeding
     differently.  Use <tt>help rand</tt> to learn a new way.
     
     <li>Possible "observations" could be on the general
     trend for 100 to 1000 to 10000 and in particular on
     the size of the
     fluctuations from a true gaussian distribution. 
     
     <li>Two of the methods are more-or-less self explanatory
     (Taylor series and diagonalizing).  The Pade approach is not.
     Be sure to look at the comments; also, a paper describing the
     approach is linked on the class webpage (at the bottom).
     You might want to google "Pade Approximant" and look at the
     Wikipedia and/or MathWorld entries.
     For part iii, the variable <tt>k</tt> corresponds to the
     number of terms.  You don't want to print it out until
     the end of the program. 
     
     
   </ol>

 <li>Fun with the One-Particle Stochastic Variational Method in MATLAB.
  
  <ol type="a">
    <li>When thinking about the advantages of a gaussian basis,
    you might consider how the matrix elements of the Hamiltonian
    are calculated here and how it would be done with a different
    basis.
    
    <li>Ill-conditioned matrices arise if the basis is close to
    being linearly dependent.  Recall that if there is linear
    dependence, the determinant will be zero.  If it is close, there
    will be at least one very small eigenvalue and so a large
    ratio of eigenvalues.  Look at <tt>cond</tt> under MATLAB help
    for a discussion of the "condition number", which is an indication
    of whether your matrix will have numerical problems from
    ill conditioning.
    
    <li>Think about the wave functions
    for the ground and first excited states in a hydrogen atom.
    How do they differ?  (Compare to the first two states for
    a square well, for example.)
    
  </ol>

 <LI>Stochastic calculation of multidimensional integrals.
 
   <OL TYPE="a">
     <li>If you get stuck on what the limits are, run the code
     a few times and look at the histogram.  How large (positive
     or negative) can x values become?
     
     <li>The difference in the histogram before and after changing
     <tt>eps</tt> to get a 50% acceptance should be significant.
     You should run the "before" case several times to appreciate
     the difference (and how much it fluctuates).
     
     <li>You can check whether thermalization is important by
     changing the value of the number of thermalization steps
     to a small number (like one).
     
     <li>What does the central limit theorem say about how the
     standard deviation scales with N?
     
     <li>The action here can be thought of as the integral over
     a potential.  What is the difference between a potential
     with positive and negative lambda?
     
   </OL>


 <LI>Continuing with the Partition Function for One Particle.
   <OL TYPE="a">
     <li>It may be useful to think of special cases for the
     matrix A.

     <li>What does the index i represent "physically"?
     
     <li>For the MATLAB part, you should think of q, q', k, and
     j as matrix indices, so the right side of Eq.(3) is just
     a matrix multiplication.  You will have to think about
     how to construct the matrices (there are short-cut ways
     using <tt>meshgrid</tt>,
     but first think about doing it with loops).
     
     <li>A check of your result is that A<sup>-1</sup>A = 1.
     
     <li>What derivatives will "bring down" <i>x<sub>i</sub></i>
     and <i>x<sub>j</sub></i> in the original Z[f] expression?

     
     
   </OL>

 <LI>Directly Solving for the Green's function G<sub>0</sub>.
   <OL TYPE="a">
    <LI>One way to proceed is to substitute for G<sup>0</sup> and
    also for the delta function and then to project out a given 
    <bf>k</bf>.  Remember that the non-interacting system is uniform.
    <LI>If you don't write these down by inspection, you're working too
    hard.  Remember that the delta function is zero in each of those
    regions.
    <LI>Try integrating the equation from tau-tau' = -epsilon to +epsilon 
    and then taking epsilon to zero.  What survives? 
   </OL>

</OL>




<HR>
<A HREF="http://www.physics.ohio-state.edu/~ntg/880_2009/">[880.05
                Home Page]</A>
<A HREF="http://www.physics.ohio-state.edu/">[OSU Physics]</A>


<HR>


<STRONG>Physics 880.05: Assignment #2 hints.</STRONG><BR>
   Last modified: <?echo Date("h:i a, F d, Y", getLastMod());?>.<BR>
  <ADDRESS><A HREF="mailto:furnstahl.1@osu.edu">furnstahl.1@osu.edu</A></ADDRESS>


</BODY>

</HTML>
