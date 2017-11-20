# Master (will become release 2.7)

# Release 2.6

*  The `diffOrder` value has disappeared from the `LocalBasisTraits` class.
   This value encoded the highest partial derivative order implemented by
   a local basis. Encoding this value as a compile-time parameter led to
   various problems related to the dynamic interface, mainly because it
   became part of the type of the local finite element.  At the same time,
   it was suspected that very few people ever actually used the parameter.

   More practically, two things have disappeared: the `diffOrder` member
   of the `LocalBasisTraits` class, and the 8th template parameter `dorder`
   of that class.  There is no replacement, and if you have used `diffOrder`
   then you currently have to find a way to live without it.  As mentioned
   we believe that this concerns only a very small number of people.

   If you do use `diffOrder` and you absolutely need it or something similar,
   then we'd like to hear from you.  One of the reasons why there is no
   replacement is that we couldn't really think of a good use case to begin with.
