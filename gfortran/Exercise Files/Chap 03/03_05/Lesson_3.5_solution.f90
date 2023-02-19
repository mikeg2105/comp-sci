program todo
    implicit none

    integer, parameter :: MAX_TASKS = 99
    integer :: i
    integer :: num_tasks
    character(len=1) :: response
    integer :: task
    character(len=100) :: tasks(MAX_TASKS)
    character(len=100) :: todo_file
    integer :: status
    integer :: unit

    select case (command_argument_count())
    case (0)
        todo_file = "todo.txt"
    case (1)
        call get_command_argument(1, todo_file)
    case default
        print *, "Too many command line arguments."
        print *, "Usage: todo [file]"
        print *, "  where file is optional, and is 'todo.txt' by default"
        stop 1
    end select

    open(newunit=unit, file=trim(todo_file))

    do i = 1, MAX_TASKS
        read(unit, '(A)', iostat=status) tasks(i)
        if (status /= 0) then
            num_tasks = i - 1
            exit
        else
            num_tasks = i
        end if
    end do

    close(unit)

    do
        if (num_tasks > 0) then
            print *, "Here are your current tasks."
            do i = 1, num_tasks
                print '(I3,") ",A)', i, trim(tasks(i))
            end do
        else
            print *, "You have no tasks."
        end if

        print *, "What would you like to do? (a)dd, (d)elete, (q)uit"
        read(*, '(A)') response

        select case (response)
        case ('a')
            if (num_tasks == MAX_TASKS) then
                print *, "Can't save any more tasks."
            else
                print *, "What's the task?"
                num_tasks = num_tasks + 1
                read(*, '(A)') tasks(num_tasks)
            end if
        case ('d')
            if (num_tasks > 0) then
                do
                    print *, "Which one would you like to delete?"
                    read(*, *, iostat = status) task
                    if (status == 0) then
                        if (task < 1) then
                            print *, "Task number must be >= 1."
                        else if (task > num_tasks) then
                            print *, "Task number must be <= ", num_tasks
                        else
                            do i = task, num_tasks-1
                                tasks(i) = tasks(i+1)
                            end do
                            num_tasks = num_tasks - 1
                            exit
                        end if
                    else
                        print *, "Sorry, I didn't understand that."
                    end if
                end do
            else
                print *, "No tasks to delete."
            end if
        case ('q')
            exit
        case default
            print *, "Sorry, I didn't understand that."
        end select

        open(newunit=unit, file=trim(todo_file), status="replace")

        do i = 1, num_tasks
            write(unit, '(A)') trim(tasks(i))
        end do

        close(unit)
    end do
end program
