Hi @lshep,

Thank you for the review! I will do my best to address these issues in my response 
below:

**RLSeq**

**vignette**


> * [ ]  Might want to show `library(dplyr)` as when I got to 4.3 I couldn't run
>   code if just following the shown coding chunks.

Thank you for pointing this out! I have now added `library(dplyr)`.

> * [ ]  Do you use parallelization in your code? Interesting after running the
>   vignette when I went to exit my R terminal I got the following message:
> 
> ```
>   > q()
> Save workspace image? [y/n/c]: n
> Error while shutting down parallel: unable to terminate some child processes
> ```
> 
> I was able to reproduce this twice.

I had not previously noticed this, but it appears that regioneR uses 
parallelization by default with > 1000 peaks. I have now updated the code to 
prevent this behavior.

**R code**
 
**general**
 
> * [ ]  getFromNamespace is documented as not to be used in production code and
>   should be avoided. What was the reasoning as its normal use case is to avoid
>   warnings about using non exported functions but the RLHub package exports
>   these functions already.

Thank you for pointing this out! I have removed all uses of `getFromNamespace`.

> * [ ]  Your using a lot of functions from packages that are not specified in the
>   NAMESPACE. Please import the package or select function appropriately in the
>   namespace from dplyr, valr, genomeinfodb, etc. If you find you are using many
>   functions from a package, the package should be in the namespace. I still
>   recommend keeping all of the `::` for completeness and to avoid conflicts of
>   functions names but having a more complete NAMESPACE is advisable.

I have now updated this! 


**featureEnrich**
 
> * [ ]  Setting seeds without users knowledge is generally discourages if there is
>   a random component that could effect reproducibility. It would be better
>   practice to have the seed as an argument to the function so a user is
>   aware. Also this entire section is concerning considering there is an issue
>   with the function in valr that needs this "workaround". Is there a different
>   package or function that could be utilized that does not have bugs?

I have now removed the use of the `bed_shuffle` function! There are no other instances 
where seeds are used.

> Please fix and respond to the above concerns. When ready please kick off a new build and request a re-review. Thank you

