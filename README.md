# SAC_Code
Attractive Surface-Based State-Attracted Control for Uncertain Nonlinear Systems

The relevant paper can be found at: W. He, Y. Yang, L. Kong, J. Peng, Y. Liu, and G. Li, “Attractive Surface-Based State-Attracted Control for Uncertain Nonlinear Systems,” IEEE Transactions on Automatic Control, doi: 10.1109/TAC.2026.3685547.

This method is mainly developed for scenarios in which PID control remains dominant in micro-scale devices and the desired control performance cannot be prescribed in advance. By introducing the concept of an attractive surface, PID control is interpreted as the superposition of two unbounded linear surfaces. Based on this perspective, an arbitrary bounded attractive-surface controller can be designed. Since the attractive surface is bounded and can be designed flexibly, the resulting tracking performance is generally superior to that of conventional PID control and can achieve prescribed control performance.

This paper presents only one possible approach for achieving prescribed performance by modifying the attractive surface. There are many possible implementation methods and functional forms for realizing this idea, so readers should not be limited to the specific form proposed in this work. The controller form presented here is mainly designed for second-order systems, with an emphasis on simplicity and ease of design. This concept can also be applied to many other control objectives, such as prescribed-time control, synchronized trajectory planning and tracking control, constraint control, and prescribed-performance-based tracking control.

If you have any theoretical or coding questions, please feel free to leave a comment.

If this paper and the developed code are helpful or inspiring to your work, please consider citing this article.


