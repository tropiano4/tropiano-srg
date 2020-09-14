<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"> 

<HTML>

<HEAD>
  <TITLE>OSU Physics: Physics 880.05 Assignment #1 Hints</TITLE>
  <LINK REV="made" HREF="mailto:webmaster@www.physics.ohio-state.edu">
  <LINK REL=STYLESHEET HREF="26x_style.css" TYPE="text/css">  
</HEAD>

<BODY BGCOLOR="#FFFFFF">


<?
   include ("physics/common_short2.php");
   physics_header("","");
?>
            

<H1>Physics 880.05: Assignment #1</H1>

Here are some hints, suggestions, and comments on the assignment.

<OL>

 <LI>MATLAB Sandbox
   <ol TYPE="a">
     <li>By "explicitly diagonalize" we
     mean calculate something like V<sup>-1</sup>AV or VAV<sup>-1</sup>
     and check if it is diagonal (check both!).   Note that the
     inverse of A is <tt>inv(A)</tt> in MATLAB (actually, this
     is just one way to calculate it).  Compare V<sup>-1</sup>
     to the Hermitian conjugate of V (given by <tt>V'</tt> in
     MATLAB). What does this tell you?
     
     <li>The word "cyclic" should appear in your answer.
     
     <li>To test whether A has to be complex or Hermitian, simply
     try out the different cases for a couple of examples each.
     
     <li>We want to know how many decimal points are good for the
     different numbers of terms.  This presumably depends on the
     matrix in detail, but you can try a few examples to see if
     there is a consistent trend.
      You may want <tt>format long</tt> to get enough digits
     in the answers.

     <li>Make sure that T and V do not commute (unlikely if they
     are chosen at random!).  You are basically looking for how
     many decimal digits improvement there is when epsilon is
     reduced by factors of 10.
     For this type of problem you should make you comparisons
     with different epsilons using the <em>same</em> T and V
     matrices so that it is clear that the differences are
     solely due to the numerical approximation.  Then you can
     try another set to see if your results are robust (or
     just the result of some accident).
   </ol>

 <li>SVM revisited.<br>
  Being stationary under arbitrary variations of the
  expansion coefficients just means that the partial derivative
  of the estimated energy with respect to any of the c<sub>i</sub>
  coefficients is equal to zero.  So just carry out this
  derivative, remembering that you have both a numerator and
  a denominator.  You'll be able to cancel one sum from the
  denominator because it is positive definite and you will
  use that the H and B matrices are symmetric.
  

 <LI>Symmetry factors.
   <OL TYPE="a">
     <li>lambda<sup>3</sup> contribution to &lt;xi<sup>2</sup>&gt;
     
     <ol type="i">
       <LI>This is just a review of the rules for diagrams (e.g., how is
       the number of vertices related to the power of lambda?).
       We're not looking for long answers!
       <LI>Follow the rules from the lecture notes.  Remember when
       calculating the vertex permutation factor (the third factor) that
       the "external lines" for &lt;xi<sup>2</sup>&gt; are nailed down.
       There are 10 distinct diagrams.
       <LI>Just plug in the numbers and add the fractions.
       <li>If you get stuck or want to check your answers,
       look at the solution 
        <a href="http://www.physics.ohio-state.edu/~ntg/880/notes/ps2ans.pdf">here
        (problem 3)</a>.
     </ol>
     
     <li>partition function with xi<sup>6</sup>
     <ol type="i">
       <li>The rule for lines should be the same as with xi<sup>4</sup>, 
       because that
       comes from the xi<sup>2</sup> term.  
       For the vertex, before we had the -lambda/4 factor times the
       combinatorical factor for four derivatives hitting j<sup>4</sup>.
       <li>How many lines should come from each vertex?
       Remember that the partition function (as opposed to its
       logarithm) includes ALL diagrams, connected and disconnected.
       The partition function has no external lines.
       <li>&lt;xi<sup>2</sup>&gt; will have two external lines.
       Only connected diagrams contribute to <li>&lt;xi<sup>2</sup>&gt;
       but what about the partition function (the question asks for Z/Z0 
       rather than ln[Z/Z0]).
       <li>Look at the solution to the previous problem 
        <a href="http://www.physics.ohio-state.edu/~ntg/880/notes/ps2ans.pdf">here
         (problem 3)</a> to see how things work.
         The exact answers via Mathematica are 1 - 5 alpha/(2 a^3)
         + 1155 alpha^2/(8 a^6) for Z/Z0 and 1/a - 15 alpha/a^4 + 1695
         alpha^2/a^7 for &lt;xi<sup>2</sup>&gt;.
       
     </ol>
     
     
     <li>replica method for &lt;xi<sup>2</sup>&gt;
     <ol type="i">
       <li>This is a bit tricky:  We are really <em>continuing</em>
       the function to n=0.  In the definition of O<sub>n</sub>
       in the problem set, the first term is the numerator we want
       for &lt;O&gt; but then we want a factor of Z/Z<sub>0</sub> in the
       <em>denominator</em>.  You can get this by multiplying 
       on top and bottom by Z/Z<sub>0</sub>.
       How many factors of Z/Z<sub>0</sub> are there now
       in the expression for O<sub>n</sub> besides the part we want?
       Write them as Z/Z<sub>0</sub> to that power.  What happens if
       you then set ("continue") n to 0?
       
       
       <li>Where do the external legs come from?  Can any index appear?
       
       <li>If you have disconnected diagrams, how many are there
       (consider what indices can appear)?
       For the connected diagrams, how many different indices are there?
       
       <li>The answer should be yes!  How would the argument change
       for another operator, such as &lt;xi<sup>4</sup>&gt;?
              
     </ol>


   </OL>


</OL>




<HR>
<A HREF="http://www.physics.ohio-state.edu/~ntg/880_2009/">[880.05
                Home Page]</A>
<A HREF="http://www.physics.ohio-state.edu/">[OSU Physics]</A>


<HR>


<STRONG>Physics 880.05: Assignment #1 hints.</STRONG><BR>
   Last modified: <?echo Date("h:i a, F d, Y", getLastMod());?>.<BR>
  <ADDRESS><A HREF="mailto:furnstahl.1@osu.edu">furnstahl.1@osu.edu</A></ADDRESS>


</BODY>

</HTML>
