class TwoThreeFourTree:
    def __init__(self, parent=None):
        self.root = []                  # create a list for a root
        self.subtrees = []              # create a list for the subtrees
        self.parent = parent            # create a parent (default is None)

    def createEmptyTree(self):
        """
        +createEmptyTree()
        This method creates an empty 2-3-4 tree.
        """
        self.root = []

    def destroy(self):
        """
        +destroy()
        Destroys a 2-3-4 tree.
        """
        del self.parent                # delete everything from the tree.
        del self.root
        del self.subtrees

    def isEmpty(self):
        """
        +isEmpty(): boolean     {query}
        This method returns True when the 2-3-4 tree is empty (no elements)
        :return: True --> when the list is empty
                 False --> when there is at least one element in the list
        """
        return len(self.root) == 0

    def givepos(self, item):
        """
        +givepos(in item: integer, out rootlenght)      {query}
        This function gives the position where the item is located in the subtrees
        :param item: The item where the position is being searched from
        :return:
        """
        if len(self.root) == 0:               # if the root is empty return 0
            return 0                          # return 0 (because there are no items in the tree)
        for i in range(len(self.root)):       # search for every value in the root
            if item <= self.root[i]:          # if the value in the root is greater than the item
                return i                      # return that value of the root
        return len(self.root)                 # if item is greater than last value of root, return de last value

    def leaf(self):
        """
        +leaf(): boolean        {query}
        This function returns True when node is a leaf and False when node is not a leaf
        """
        return len(self.subtrees) == 0

    def split(self):
        """
        +split(in TwoThreeFourTree: class, out TwoThreeFourTree: class)
        This function splits a node on the path of the item only if the node on that path is full (3 items)
        :return: Returns a class TwoThreeFourTree
        """
        if self.parent is None:                                         # if split for the root
            lefttree, righttree = TwoThreeFourTree(), TwoThreeFourTree()  # make two subtrees -> lefttree and righttree
            lefttree.parent, righttree.parent = self, self              # parent of subtrees is the root of self
            if not self.leaf():                                         # if the root is not a leaf
                lefttree.subtrees.append(self.subtrees[0])
                lefttree.subtrees.append(self.subtrees[1])
                self.subtrees[0].parent = lefttree                      # parent of subtrees by value is lefttree
                self.subtrees[1].parent = lefttree
                righttree.subtrees.append(self.subtrees[2])
                righttree.subtrees.append(self.subtrees[3])
                self.subtrees[2].parent = righttree                     # parent of subtrees is value of righttree
                self.subtrees[3].parent = righttree
            lefttree.root.append(self.root[0])                          # set smallest value of the root in the lefttree
            righttree.root.append(self.root[2])                         # set largest value of the root in the righttree
            self.subtrees.clear()                                       # remove all the elements in the subtrees
            self.root.pop(2)                                            # remove the largest value of the root
            self.root.pop(0)                                            # remove the smallest value of the root
            self.subtrees.append(lefttree)                              # append the lefttree to the subtrees
            self.subtrees.append(righttree)                             # append the righttree to the subtrees
        else:                                                           # if not a root (intern node or leaf)
            position = self.parent.givepos(self.root[1])
            self.parent.root.insert(position, self.root[1])
            tree = TwoThreeFourTree()
            tree.root.append(self.root[2])
            tree.parent = self.parent
            self.root.pop(2)
            self.root.pop(1)
            if not self.leaf():
                tree.subtrees.append(self.subtrees[2])
                tree.subtrees.append(self.subtrees[3])
                tree.subtrees[0].parent = tree
                tree.subtrees[1].parent = tree
                self.subtrees.pop(3)
                self.subtrees.pop(2)
            self.parent.subtrees.insert(position+1, tree)
        if self.parent:
            return self.parent
        return self

    def insert(self, item):
        """
        +insert(in TwoThreeFourTree: self, in item: class/list/integer, out success: boolean)
        Adds new item to a 2-3-4 tree with items with different
        search key values, different from the search key of new item.
        Success indicates whether the addition was successful.
        :param item: The item that's going to be inserted, different from the items that are already in the tree.
        :return: Returns the item that is inserted and success = True if item is succesfully inserted. Success = False
        if item is not succesfully inserted.
        """
        if not self.root:                                               # if root is empty
            self.root.append(item)                                      # append item to the root
            return item, True
        else:                                                           # if root is not empty
            current = self                                              # current is self
            while not current.leaf():                                   # while current is not a leaf
                if len(current.root) == 3:                              # if current's root is full
                    current.split()                                     # split this root
                    if current.parent is not None:
                        current = current.parent
                current = current.subtrees[current.givepos(item)]       # set current to the current his right subtree

            if len(current.root) != 3:                                  # if current is leaf and current's root not full
                current.root.insert(current.givepos(item), item)        # insert the item in the current root
                return item, True
            else:                                                       # else if current's root is full
                current = current.split()                               # split current, this is the new current
                current = current.subtrees[current.givepos(item)]
                current.root.insert(current.givepos(item), item)
                return item, True

    def retrieve(self, item):
        """
        +retrieve(in TwoThreeFourTree: class, in item: class/list/integer, out success: boolean)
        :param item: The item that is going to be retrieved.
        :return: A boolean: True if the item is retrieved in the Binary Tree.
                            False if the item is not retrieved in the Binary Tree.
        """
        if item not in self.root:
            if len(self.subtrees) == 0:
                return item, False
            if len(self.subtrees) == 2:
                if item < self.root[0]:
                    return self.subtrees[0].retrieve(item)
                else:
                    return self.subtrees[1].retrieve(item)
            if len(self.subtrees) == 3:
                if item < self.root[0]:
                    return self.subtrees[0].retrieve(item)
                elif item < self.root[1]:
                    return self.subtrees[1].retrieve(item)
                else:
                    return self.subtrees[2].retrieve(item)
            if len(self.subtrees) == 4:
                if item < self.root[0]:
                    return self.subtrees[0].retrieve(item)
                elif item < self.root[1]:
                    return self.subtrees[1].retrieve(item)
                elif item < self.root[2]:
                    return self.subtrees[2].retrieve(item)
                else:
                    return self.subtrees[3].retrieve(item)
        else:
            return item, True

    def merge(self):
        """
        +merge(inout TwoThreeFourTree: class)
        This function merges nodes together on the path of the item that is going to be deleted
        """
        position = self.parent.givepos(self.root[0])
        if position == 0:
            self.root.append(self.parent.root[position])
            self.root.append(self.parent.subtrees[position+1].root[0])
            if not self.leaf():
                self.subtrees.append(self.parent.subtrees[position+1].subtrees[0])
                self.subtrees.append(self.parent.subtrees[position+1].subtrees[1])
            self.parent.subtrees.pop(position+1)
            self.parent.root.pop(0)
            return
        self.root.insert(0, self.parent.root[position-1])
        self.root.insert(0, self.parent.subtrees[position-1].root[0])
        if not self.parent.subtrees[position-1].leaf():
            self.subtrees.insert(0, self.parent.subtrees[position-1].subtrees[1])
            self.subtrees.insert(0, self.parent.subtrees[position-1].subtrees[0])
        self.parent.subtrees.pop(position-1)
        self.parent.root.pop(position-1)

    def redistribute(self):
        """
        +redistribute(in TwoThreeFourTree: class, out item: integer): boolean
        This function redistributes the item on a path for deleting aan item.
        :return: An item (integer), this is the item that is deleted.
                 A boolean: True if item is deleted out of the TwoThreeFourTree
                            False if item is not deleted out of the TwoThreeFourTree
        """
        if self.parent is None:
            return True
        position = self.parent.givepos(self.root[0])
        if position == 0:
            if len(self.parent.subtrees[position+1].root) != 1:
                if not self.parent.subtrees[position+1].leaf():
                    self.subtrees.append(self.parent.subtrees[position+1].subtrees[0])
                    self.parent.subtrees[position+1].subtrees.pop(0)
                self.root.append(self.parent.root[0])
                self.parent.root[0] = self.parent.subtrees[position+1].root[0]
                self.parent.subtrees[position+1].root.pop(0)
                return True

        elif position == len(self.parent.root):
            if len(self.parent.subtrees[position - 1].root) != 1:
                if not self.parent.subtrees[position - 1].leaf():
                    self.subtrees.insert(0, self.parent.subtrees[position - 1].subtrees[
                        len(self.parent.subtrees[position - 1].subtrees) - 1])
                    self.parent.subtrees[position - 1].subtrees.pop(
                        len(self.parent.subtrees[position - 1].subtrees) - 1)
                self.root.insert(0, self.parent.root[len(self.parent.root) - 1])
                self.parent.root[len(self.parent.root) - 1] = self.parent.subtrees[position - 1].root[
                    len(self.parent.subtrees[position - 1].root) - 1]
                self.parent.subtrees[position - 1].root.pop(-1)
                return True
            else:
                return False
        else:
            if len(self.parent.subtrees[position - 1].root) == 1:
                if len(self.parent.subtrees[position + 1].root) != 1:
                    if not self.parent.subtrees[position + 1].leaf():
                        self.subtrees.append(self.parent.subtrees[position + 1].subtrees[0])
                        self.parent.subtrees[position + 1].subtrees.pop(0)
                    self.root.append(self.parent.root[position])
                    self.parent.root[position] = self.parent.subtrees[position + 1].root[0]
                    self.parent.subtrees[position + 1].root.pop(0)
                    return True
                else:
                    return False
            else:
                if not self.parent.subtrees[position-1].leaf():
                    self.subtrees.insert(0, self.parent.subtrees[position-1].subtrees[
                        len(self.parent.subtrees[position-1].subtrees)-1])
                    self.parent.subtrees[position-1].subtrees.pop(len(self.parent.subtrees[position-1].subtrees)-1)
                self.root.insert(0, self.parent.root[position - 1])
                self.root.sort()
                self.parent.root[position-1] = self.parent.subtrees[position-1].root[-1]
                self.parent.subtrees[position-1].root.pop(-1)
                return True

    def inordersuccesor(self, item):
        """
        +inordersuccesor(in item: class/list/integer, out current: class)      {query}
        This function gives the inorder succesor of an item in a node.
        :param item: This parameter is an integer and this is the item where the inorder succesor has to be found.
        :return: Returns the node where the inorder succesor is found.
        """
        current = self.subtrees[self.givepos(item)+1]       # go to the lefttree of the item
        if len(current.root) == 1:
            print("NOG DOEN")
        while not current.leaf():                           # while the current is not a leaf
            current = current.subtrees[0]                   # go to the most right tree
            if not current.leaf() and len(current.root) == 1:
                print("NOG DOEN")
        return current                                      # if current is a leaf than return current

    def inorderTraversal(self, array=None):
        """
        +inorderTraversal(in TwoThreeFourTree: class, in array: vector of integers)
        This function traverses every item in the tree in the order: first go to the lefttree of the TwoThreeFourTree.
        If lefttree fully traversed than go to the root and after the root, go to the righttree of the TwoThreeFourTree.
        :param array: array is an empty list. This list is the list to keep the order of the items in which they are
        travesed.
        :return: The list "lijst" (see in param)
        """
        if array is None:                   # make an empty array, this array holds the values who are sorted
            array = []
        for i in self.root:                 # loop over every value in the root and append this to the array
            array.append(i)
        for j in self.subtrees:             # loop over every subtree and apply inorderTraversal
            j.inorderTraversal(array)
        array.sort()
        return array                        # return array

    def delete(self, item):
        """
        +delete(in item: class/list/integer, out class/list/integer, out success: boolean)
        This function deletes an item in a 2-3-4 Tree and returns the item and a boolean.
        :param item: This is the item that is going to be deleted in the 2-3-4 Tree.
        :return: Returns the item that is deleted and a boolean: True if item is deleted
                                                                 False if item is not deleted
        """
        if not self.retrieve(item)[1]:                                                  # Als item niet in boom zit
            return item, False
        else:
            current = self                                                              # Zet current op self
            if len(self.root) == 1:                                                     # Als de root één item heeft
                if len(self.subtrees[0].root) == 1 and len(self.subtrees[1].root) == 1: # Beide subtrees één item hebben
                    self.root.insert(0, self.subtrees[0].root[0])                       # Dit is de merge voor de root
                    self.root.append(self.subtrees[1].root[0])
                    if not self.subtrees[0].leaf():
                        self.subtrees.append(self.subtrees[0].subtrees[0])
                        self.subtrees.append(self.subtrees[0].subtrees[1])
                        self.subtrees.append(self.subtrees[1].subtrees[0])
                        self.subtrees.append(self.subtrees[1].subtrees[1])
                    self.subtrees.pop(0)
                    self.subtrees.pop(0)
                    for i in self.subtrees:
                        i.parent = self

            while item not in current.root:                                            # while item not in current root
                if len(current.subtrees[current.givepos(item)].root) == 1:             # Als root van subtree van item 1
                    if not current.subtrees[current.givepos(item)].redistribute():     # Als redistribute niet gaat
                        current.subtrees[current.givepos(item)].merge()                # Doe dan een merge
                current = current.subtrees[current.givepos(item)]                      # verander current tot root item

            if current.leaf():                                                         # if current is a leaf
                current.root.remove(item)                                              # remove item out of leaf
                return item, True

            else:                                                                       # Dit is voor een interne node of root als de inorder succ moet worden gevonden
                hold = current.givepos(item)
                if len(current.subtrees[hold].root) > 1:
                    current.root[hold] = current.subtrees[hold].root[-1]
                    del current.subtrees[hold].root[-1]
                else:
                    suc = current.inordersuccesor(item)
                    hold = current.givepos(item)
                    current.root[current.givepos(item)] = suc.root[0]
                    thesucitem = suc.root[0]
                    if len(suc.root) == 1:
                        if len(suc.parent.root) == 1:
                            suc.root.pop()
                            suc.root.append(suc.parent.root[0])
                            del suc.parent.root[0]
                            suc.parent.root.append(suc.parent.subtrees[1].root[0])
                            del suc.parent.subtrees[1].root[0]
                        elif len(suc.parent.root) == 2:
                            if hold == 0:
                                suc.parent.subtrees[hold].merge()
                                suc.parent.subtrees[hold].root.remove(thesucitem)
                            elif hold == 1:
                                current.root.pop(1)
                                current.subtrees[1].root.extend(current.subtrees[2].root)
                                del current.subtrees[2]
                        elif len(suc.parent.root) == 3:
                            if hold == 0:
                                suc.parent.subtrees[0].merge()
                                suc.parent.subtrees[0].root.remove(thesucitem)
                            elif hold == 1:
                                suc.parent.subtrees[0].merge()
                                suc.parent.subtrees[0].root.remove(thesucitem)
                            elif hold == 2:
                                current.root.pop(2)
                                current.subtrees[2].root.extend(current.subtrees[3].root)
                                del current.subtrees[3]
                    elif len(suc.root) == 2:
                        suc.root.remove(thesucitem)

    def content(self, file):
        """
        +content(in file: .dot-file)    {query}
        This function is an auxiliary function for the function "toDot". This function writes the content
        of the .dot file.
        """
        if len(self.subtrees) == 0:
            if len(self.root) == 1:
                file.writelines("\"" + str(self.root[0]) + "\"")
            elif len(self.root) == 2:
                file.writelines("\"" + str(self.root[0]) + " | " + str(self.root[1]) + "\"")
            else:
                file.writelines("\"" + str(self.root[0]) + " | " + str(self.root[1]) + " | " + str(self.root[2]) + "\"")
        else:
            for i in self.subtrees:
                file.write("\"")
                if len(self.root) == 1:
                    file.writelines(str(self.root[0]))
                elif len(self.root) == 2:
                    file.writelines(str(self.root[0]) + " | " + str(self.root[1]))
                else:
                    file.writelines(str(self.root[0]) + " | " + str(self.root[1]) + " | " + str(self.root[2]))
                file.write("\"--\"")
                if len(i.root) == 1:
                    file.writelines(str(i.root[0]))
                elif len(i.root) == 2:
                    file.writelines(str(i.root[0]) + " | " + str(i.root[1]))
                else:
                    file.writelines(str(i.root[0]) + " | " + str(i.root[1]) + " | " + str(i.root[2]))
                file.write("\" \n")
                i.content(file)

    def toDot(self, name):
        """
        +toDot(in name: string)    {query}
        This function writes a .dot file with a given name for the TwoThreeFourTree.
        """
        file = open(str(name) + ".dot", "w+")       # make a new (or overwrite existing).dot-file that has a name
        file.write("Graph{\n")                      # first line always begins with "Graph{" + newline
        self.content(file)                          # write the contents of the file (this is every node of the tree)
        file.write("}")                             # end file always with "}"


a = TwoThreeFourTree()
a.insert(1)
a.insert(2)
a.insert(3)
a.insert(4)
#a.delete(4)
#a.delete(1)
a.insert(5)
a.insert(6)
a.insert(7)
a.insert(8)
a.insert(9)
a.insert(10)
a.insert(11)
#a.delete(3)
a.insert(12)
a.insert(13)
a.insert(14)
a.insert(15)
a.insert(16)
a.insert(17)
a.insert(18)
a.insert(19)
a.insert(20)
#a.delete(8)
#a.delete(8)
#a.delete(9)
#a.delete(10)
#a.delete(14)
#a.delete(12)
#a.delete(6)
#a.delete(3)
#a.delete(11)
#a.delete(13)
#a.delete(12)
#a.insert(15)
print(a.inorderTraversal())

###################
#a.insert(1)
#a.insert(10)
#a.insert(3)
#a.insert(8)
#a.insert(5)
#a.delete(1)
#a.insert(6)
#a.insert(7)
#a.insert(4)
#a.delete(8)
#a.insert(9)
#a.insert(2)
#a.delete(3)
#a.delete(7)
a.toDot("EY")